library(magrittr)
library(tidyr)
library(dplyr)

observed_data_loglik <- function(n, n_v, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data_mat, theta, p, gamma, q, NP_YstarDist)
{
  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[c(1:n_v),c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(data_mat[c(1:n_v),c(Y_val)]) == 0, 1-pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over log[P(Yi*|Xi*,Yi,Xi)] -----------------------------------------------
  if (NP_YstarDist)
  {
    suppressMessages(phat <- data_mat[c(1:n_v),] %>% data.frame(row.names = NULL) %>% 
      dplyr::left_join(data.frame(p)) %>% 
      dplyr::mutate(pYstar = ifelse(is.na(pYstar), 0, pYstar)) %>% 
      dplyr::pull(pYstar))
    log_phat <- ifelse(phat == 0, 0, log(phat))
    return_loglik <- return_loglik + sum(log_phat)
  } else
  {
    phat <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[c(1:n_v),c(X_unval, Y_val, X_val, C)]) %*% gamma)))
    phat <- ifelse(as.vector(data_mat[c(1:n_v),Y_unval]) == 0, 1-phat, phat)
    return_loglik <- return_loglik + sum(log(phat))
  }
  ## ----------------------------------------------- Sum over log[P(Yi*|Xi*,Yi,Xi)]
  #################################################################################
  ## Sum over I(Xi=xk)Bj(Xi*)log q_kj ---------------------------------------------
  sn <- ncol(q)
  m <- nrow(q)
  log_q <- log(q)
  log_q[log_q == -Inf] <- 0
  unique_X <- data_mat[1:n_v,] %>% data.frame() %>% dplyr::select(X_val) %>% unique()
  log_q %<>% data.frame() %>% magrittr::set_colnames(paste0("log_q_bs", seq(1,sn))) %>% 
    cbind(unique_X)
  suppressMessages(val_data_mat <- data_mat[1:n_v,] %>% data.frame() %>% dplyr::left_join(log_q) %>% data.matrix())
  return_loglik <- return_loglik + sum(val_data_mat[,Bspline] * val_data_mat[,paste0("log_q_bs", seq(1,sn))])
  ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  #################################################################################
  # -------------------------------------------------------- For validated subjects
  
  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[-c(1:n_v),c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(data_mat[-c(1:n_v),c(Y_val)]) == 0, 1-pY_X, pY_X)
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate P(Y*|X*,y,x) for all (y,xk) ---------------------------------------
  if (NP_YstarDist)
  {
    suppressMessages(phat <- data_mat[-c(1:n_v),] %>% data.frame(row.names = NULL) %>% 
                       dplyr::left_join(data.frame(p)) %>% 
                       dplyr::mutate(pYstar = ifelse(is.na(pYstar), 0, pYstar)) %>% 
                       dplyr::pull(pYstar))
  } else
  {
    phat <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[-c(1:n_v),c(X_unval, Y_val, X_val, C)]) %*% gamma)))
    phat <- ifelse(as.vector(data_mat[-c(1:n_v),Y_unval]) == 0, 1-phat, phat)
  }
  ## --------------------------------------- Calculate P(Y*|X*,y,x) for all (y,xk)
  ################################################################################
  ## Merge in q_kj ---------------------------------------------------------------
  q %<>% data.frame() %>% magrittr::set_colnames(paste0("q_bs", seq(1,sn))) %>% 
    cbind(unique_X)
  suppressMessages(unval_data_mat <- data_mat[-c(1:n_v),] %>% data.frame() %>% dplyr::left_join(q) %>% data.matrix())
  ## --------------------------------------------------------------- Merge in q_kj
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x q_kj ---------------------------
  inside_product <- c(pY_X*ifelse(is.na(phat), 0, phat))*unval_data_mat[,Bspline] * unval_data_mat[,paste0("q_bs", seq(1,sn))]
  inside_product_sum <- rowSums(rowsum(inside_product, group = rep(seq((n_v+1), n), times = 2*m)))
  ## And sum over log of them all ------------------------------------------------
  return_loglik <- return_loglik + sum(ifelse(inside_product_sum == 0, 0, log(inside_product_sum)))
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}

profile_out <- function(theta, n_v, n, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data_mat, p0, gamma0, q0, q_val_num, q_val_denom, NP_YstarDist, VERBOSE = FALSE, TOL = 1E-4, MAX_ITER = 5000)
{
  sn <- ncol(q0)
  m <- nrow(q0)
  
  prev_pYstar <- p0
  prev_gamma <- gamma0
  prev_q <- q0
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- ""
  it <- 1
  
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[,c(X_val, C)]) %*% theta)))
    pY_X <- ifelse(as.vector(data_mat[,c(Y_val)]) == 0, 1-pY_X, pY_X)
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    if (NP_YstarDist)
    {
      suppressMessages(phat <- data_mat %>% data.frame(row.names = NULL) %>% 
                         dplyr::left_join(data.frame(prev_pYstar, row.names = NULL)) %>% dplyr::pull(pYstar))
      phat[is.na(phat)] <- 0
    } else
    {
      phat <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)]) %*% prev_gamma)))
      phat <- ifelse(as.vector(data_mat[,Y_unval]) == 0, 1-phat, phat)
    }
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### q_kj ----------------------------------------------------------
    prev_q_stacked <- matrix(0, nrow = n_v, ncol = sn) %>% 
      rbind(prev_q[rep(seq(1,nrow(prev_q)), each = 2*(n-n_v)),])
    ### ---------------------------------------------------------- q_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)B(X*)log(q) ----------------------------
    psi_mat_t_num <- pY_X * phat * data_mat[,Bspline] * prev_q_stacked
    ### Update denominator --------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ------------------
    psi_mat_t_rowsums <- rowsum(psi_mat_t_num, group = c(seq(1,n_v), rep(seq((n_v+1), n), times = 2*m)))
    #### Then sum over the sn splines ---------------------------------
    psi_mat_t_denom <- rowSums(psi_mat_t_rowsums)
    psi_mat_t_denom_stacked <- matrix(rep(c(rep(1,n_v), rep(psi_mat_t_denom[-c(1:n_v)], times = (2*m))), times = sn), ncol =sn, byrow = FALSE)
    ### And divide them! ----------------------------------------------
    psi_mat_t <- psi_mat_t_num/psi_mat_t_denom_stacked
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
    ## Update the Psi_kyi for unvalidated subjects --------------------
    ## by summing across the splines/ columns of psi_mat_t ------------
    Psi_mat_t <- rowSums(psi_mat_t)
    ## For validated subjects, Psi_kyi = I(Xi=xk, Yi=y) ---------------
    Psi_mat_t[1:n_v] <- 1
    ## -------------------- Update the Psi_kyi for unvalidated subjects
    ## Check that they sum to 1 
    ## data.frame(id = c(seq(1,n_v), rep(seq((n_v+1), n), times = 2*m)), Psi_mat_t) %>% group_by(id) %>% summarize(sumPsi = sum(Psi_mat_t)) %>% arrange(desc(sumPsi))
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update P(Y*|X*,Y,X,C) ------------------------------------------
    if (NP_YstarDist)
    {
      new_pYstar <- p0
      for (i in 1:nrow(new_pYstar))
      {
        which_Yunval <- as.vector(which(data_mat[-c(1:n_v),Y_unval] == new_pYstar[i,Y_unval]))
        which_Yval <- as.vector(which(data_mat[-c(1:n_v),Y_val] == new_pYstar[i,Y_val]))
        
        # allow for multiple Xs 
        which_Xunval <- as.vector(which(data_mat[-c(1:n_v),X_unval[1]] == new_pYstar[i,X_unval[1]]))
        which_Xval <- as.vector(which(data_mat[-c(1:n_v),X_val[1]] == new_pYstar[i,X_val[1]]))
        if (length(X_val) > 1)
        {
          for (d in 2:length(X_val))
          {
            which_Xunval %<>% intersect(as.vector(which(data_mat[-c(1:n_v),X_unval[d]] == new_pYstar[i,X_unval[d]])))
            which_Xval %<>% as.vector(which(data_mat[-c(1:n_v),X_val[d]] == new_pYstar[i,X_val[d]]))
          }
        }
        new_pYstar[i,"Num_Freq"] <- new_pYstar[i,"Num_Freq"] + sum(Psi_mat_t[intersect(intersect(intersect(which_Xunval, which_Yval), which_Xval), which_Yunval)])
        new_pYstar[i,"Denom_Freq"] <- new_pYstar[i,"Denom_Freq"] + sum(Psi_mat_t[intersect(intersect(which_Xunval, which_Yval), which_Xval)])
      }
      new_pYstar[,"pYstar"] <- new_pYstar[,"Num_Freq"]/new_pYstar[,"Denom_Freq"]
      #if(VERBOSE) print(new_p)
      # Check for convergence -----------------------------------------
      p_conv <- abs(new_pYstar[,"pYstar"] - prev_pYstar[,"pYstar"]) < TOL
      if(VERBOSE & mean(p_conv) == 1) print("{P(Y*|X*,Y,X) coverged")
      new_gamma <- NULL
      gamma_conv <- 1
    } else
    {
      ### Gradient ----------------------------------------------------
      gamma_design_mat <- cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)])
      mu <- gamma_design_mat %*% prev_gamma
      gradient_gamma <- matrix(data = c(colSums(Psi_mat_t * c((data_mat[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
      ### Hessian -----------------------------------------------------
      hessian_gamma <- matrix(0, nrow = ncol(gamma_design_mat), ncol = ncol(gamma_design_mat), byrow = TRUE)
      for (l in 1:ncol(gamma_design_mat))
      {
        hessian_gamma[l,] <- colSums(c(Psi_mat_t*gamma_design_mat[,l])*gamma_design_mat*c((exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1)))
      }
      new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma, 
                            error = function(err) {matrix(NA, nrow = nrow(prev_gamma))})
      if (TRUE %in% is.na(new_gamma))
      {
        warning("Hessian(gamma) is not invertible. Iteration stopped before convergence.")
        CONVERGED <- FALSE
        CONVERGED_MSG <- "Hessian gamma is not invertible"
        return(list("psi_mat_at_conv" = psi_mat_t,
                    "p_at_conv" = new_pYstar, 
                    "gamma_at_conv" = new_gamma,
                    "q_at_conv" = new_q,
                    converged = CONVERGED,
                    converged_message = CONVERGED_MSG))
      }
      p_conv <- 1
      new_pYstar <- NULL
      gamma_conv <- abs(new_gamma - prev_gamma)<TOL
      if(VERBOSE & mean(gamma_conv) == 1) print("Gamma converged")
      if(VERBOSE) print(new_gamma)
    }
    ## ------------------------------------------ Update P(Y*|X*,Y,X,C)
    ###################################################################
    ## Update {q_kj} --------------------------------------------------
    ### Update numerators ---------------------------------------------
    new_q_num <- q_val_num + rowsum(psi_mat_t[-c(1:n_v),], group = rep(seq(1,m), each = 2*(n-n_v)))
    ### Update denominators -------------------------------------------
    new_q_denom <- q_val_denom + matrix(data = rep(colSums(psi_mat_t[-c(1:n_v),]), each = m), nrow = m, ncol = sn)
    ### And finally divide them! --------------------------------------
    new_q <- new_q_num/new_q_denom
    #if(VERBOSE) print(new_q)
    ### Check for convergence -----------------------------------------
    #q_conv <- abs(new_q - prev_q)<TOL
    q_conv <- abs(new_q - prev_q)<1E-3
    ## -------------------------------------------------- Update {q_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(VERBOSE) print(paste("Iteration", it, "complete."))
    
    all_conv <- c(p_conv, gamma_conv, q_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_pYstar <- new_pYstar
    prev_gamma <- new_gamma
    prev_q <- new_q
  }
  if (!CONVERGED & it >= MAX_ITER)
  {
    CONVERGED_MSG <- "MAX_ITER reached."
  }
  return(list("psi_mat_at_conv" = psi_mat_t,
              "p_at_conv" = new_pYstar, 
              "gamma_at_conv" = new_gamma,
              "q_at_conv" = new_q,
              converged = CONVERGED,
              converged_message = CONVERGED_MSG))
}

TwoPhase_LogReg <- function(Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, Val = NULL, data, h_n_scale = 1, noSE = FALSE, VERBOSE = FALSE, TOL = 1E-4, MAX_ITER = 1000, NP_YstarDist = NULL, NP_threshold = 0)
{
  n <- nrow(data)
  V <- data[,Val] # validation indicators
  n_v <- sum(V)
  id <- seq(1,n)
  
  # Reorder so that the n_v validated subjects are first ------------
  data <- data[order(as.numeric(V), decreasing = TRUE),]
  
  # If unvalidated variable was left blank, assume error-free -------
  if (is.null(Y_unval)) Y_unval <- Y_val
  if (is.null(X_unval)) X_unval <- X_val
  
  Ystar <- data %>% dplyr::select(Y_unval)
  Y <- data %>% dplyr::select(Y_val)
  Xstar <- data %>% dplyr::select(X_unval)
  X <- data %>% dplyr::select(X_val)
  if(!is.null(C)) Covar <- data %>% dplyr::select(C)
  
  # Add the B spline basis ------------------------------------------
  B <- data[,Bspline]
  sn <- ncol(B)
  if(0 %in% colSums(B[1:n_v,])) 
  {
    warning("The B-spline basis provided has empty sieves in the Phase II subsample. Please re-specify and try again.")
    return(list(Coefficients = data.frame(Coefficient = rep(NA, length(X_val)+1), 
                                          SE = NA),
                converged = FALSE,
                converged_message = "Empty sieve in B-spline basis"))
  }
  # ------------------------------------------ Add the B spline basis
  
  # Add m for distinct values of X ----------------------------------
  x_k <- data[1:n_v,] %>% dplyr::select(X_val) %>% unique()
  m <- nrow(x_k)
  # ---------------------------------- Add m for distinct values of X 
  
  p <- length(X_val) + length(C) # number of coefficients in the analysis model (beta_1, ..., beta_p)
  prev_theta <- matrix(0, nrow = (p+1), ncol = 1)
  
  # If a modeling preference for P(Y*|X*,Y,X) -----------------------
  # wasn't supplied, test for it ------------------------------------
  if (is.null(NP_YstarDist))
  {
    NP_YstarDist <- ifelse(mean(data.frame(table(data[c(1:n_v), c(X_unval, Y_val, X_val, C)]))$Freq == 0) > NP_threshold, 
                           FALSE, 
                           TRUE)
  }
  
  if(NP_YstarDist)
  {
    tab_all <- data[c(1:n_v),c(Y_unval, X_unval, Y_val, X_val, C)] %>% dplyr::group_by_all() %>% count() %>% data.frame()
    colnames(tab_all)[length(colnames(tab_all))] <- "Num_Freq"
    tab_exclYstar <- data[c(1:n_v),c(X_unval, Y_val, X_val, C)] %>% dplyr::group_by_all() %>% count() %>% data.frame()
    colnames(tab_exclYstar)[length(colnames(tab_exclYstar))] <- "Denom_Freq"
    suppressMessages(p0 <- prev_pYstar <- tab_all %>% dplyr::full_join(tab_exclYstar) %>% 
                       dplyr::mutate(pYstar = Num_Freq/Denom_Freq) %>% data.matrix()) 
    gamma0 <- prev_gamma <- NULL
  } else
  {
    #warning("Logistic regression will be used to estimate P(Y*|X*,Y,X,C).")
    NP_YstarDist <- FALSE
    prev_gamma <- gamma0 <- matrix(0, nrow = (length(X_unval) + length(Y_val) + length(X_val) + length(C) + 1), ncol = 1)
    # cd_gamma_df <- cbind(Ystar, Xstar, Y, X) 
    # if(!is.null(C)) cd_gamma_df %<>% cbind(Covar)
    # colnames(cd_gamma_df)[1] <- "Ystar"
    # colnames(cd_gamma_df)[-1] <- paste0("X", seq(1,ncol(cd_gamma_df)-1))
    # cd_gamma <- glm(Ystar ~ ., family = "binomial", data = cd_gamma_df)
    # prev_gamma <- gamma0 <- matrix(cd_gamma$coefficients, ncol = 1)
    p0 <- prev_pYstar <- NULL
  }

  # Replace initial q w/ Bspline version  ------------ 
  prev_q <- q_val_num <- q_val_denom <- matrix(0, nrow = m, ncol = sn)
  for (j in 1:sn)
  {
    for (k in 1:m)
    {
      which_i <- which(X[1:n_v,1] == matrix(rep(x_k[k,1], n_v), nrow = n_v, byrow = TRUE))
      if (length(X_val) > 1)
      {
        for (d in 2:length(X_val))
        {
          which_i %<>% intersect(which(X[1:n_v,d] == matrix(rep(x_k[k,d], n_v), nrow = n_v, byrow = TRUE)))
        }
      }
      prev_q[k,j] <- q_val_num[k,j] <- sum(B[which_i,j])
    }
    q_val_denom[,j] <- sum(prev_q[,j]) 
    prev_q[,j] <- prev_q[,j]/sum(prev_q[,j])
  }
  colnames(prev_q) <- paste0("q_bs", seq(1,sn))
  q0 <- prev_q
  # check that q_kj sum to 1 over each j (column/ spline)
  # TRUE %in% round(colSums(q0),6) != 1
  # ------------ Replace initial q w/ Bspline version 
  
  # Add the Bspline to the psi_mat_t ----------------
  ## For validated subjects, keep only Xi/Yi --------
  ## (So one row per person) ------------------------
  ## For unvalidated subjects, create a row for -----
  ## each observed (xk,y) ---------------------------
  if(!is.null(C))
  {
    # validated subjects   
    data_mat <- data.matrix(cbind(X, Y, Xstar, Ystar, Covar, B)[1:n_v,])
    # unvalidated subjects (all combinations we might have observed)
    x_rep <- matrix(rep(x_k[,1], each = 2*(n-n_v)), ncol = 1, byrow = TRUE)
    if (ncol(x_k) > 1)
    {
      for (d in 2:ncol(x_k))
      {
        x_rep %<>% cbind(matrix(rep(x_k[,d], each = 2*(n-n_v)), ncol = 1, byrow = TRUE))
      }
    }
    colnames(x_rep) <- X_val
    suppressWarnings(data_mat %<>% rbind(data.matrix(cbind(x_rep, 
                                                           Y = rep(rep(c(0,1), each = (n-n_v)), times = m),
                                                           data[-c(1:n_v),c(X_unval, Y_unval, C)], 
                                                           B[-c(1:n_v),]))))
  } else 
  {
    # validated subjects   
    data_mat <- data.matrix(cbind(X, Y, Xstar, Ystar, B)[1:n_v,])
    # unvalidated subjects (all combinations we might have observed)
    x_rep <- matrix(rep(x_k[,1], each = 2*(n-n_v)), ncol = 1, byrow = TRUE)
    if (ncol(x_k) > 1)
    {
      for (d in 2:ncol(x_k))
      {
        x_rep %<>% cbind(matrix(rep(x_k[,d], each = 2*(n-n_v)), ncol = 1, byrow = TRUE))
      }
    }
    colnames(x_rep) <- X_val
    suppressWarnings(data_mat %<>% rbind(data.matrix(cbind(x_rep, 
                                          Y = rep(rep(c(0,1), each = (n-n_v)), times = m),
                                          data[-c(1:n_v),c(X_unval, Y_unval)], 
                                          B[-c(1:n_v),]))))
  }
  # ---------------- Add the Bspline to the psi_mat_t
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- ""
  it <- 1

  # Estimate theta using EM ----------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[,c(X_val, C)]) %*% prev_theta)))
    pY_X <- ifelse(as.vector(data_mat[,c(Y_val)]) == 0, 1-pY_X, pY_X)
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    if (NP_YstarDist)
    {
      suppressMessages(phat <- data_mat %>% data.frame(row.names = NULL) %>% 
                         dplyr::left_join(data.frame(prev_pYstar, row.names = NULL)) %>% dplyr::pull(pYstar))
      phat[is.na(phat)] <- 0
    } else
    {
      phat <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)]) %*% prev_gamma)))
      phat <- ifelse(as.vector(data_mat[,Y_unval]) == 0, 1-phat, phat)
    }
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### q_kj ----------------------------------------------------------
    prev_q_stacked <- matrix(0, nrow = n_v, ncol = sn) %>% 
      rbind(prev_q[rep(seq(1,nrow(prev_q)), each = 2*(n-n_v)),])
    ### ---------------------------------------------------------- q_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)B(X*)log(q) ----------------------------
    psi_mat_t_num <- pY_X * phat * data_mat[,Bspline] * prev_q_stacked
    ### Update denominator --------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ------------------
    psi_mat_t_rowsums <- rowsum(psi_mat_t_num, group = c(seq(1,n_v), rep(seq((n_v+1), n), times = 2*m)))
    #### Then sum over the sn splines ---------------------------------
    psi_mat_t_denom <- rowSums(psi_mat_t_rowsums)
    psi_mat_t_denom_stacked <- matrix(rep(c(rep(1,n_v), rep(psi_mat_t_denom[-c(1:n_v)], times = (2*m))), times = sn), ncol =sn, byrow = FALSE)
    ### And divide them! ----------------------------------------------
    psi_mat_t <- psi_mat_t_num/psi_mat_t_denom_stacked
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
    ## Update the Psi_kyi for unvalidated subjects --------------------
    ## by summing across the splines/ columns of psi_mat_t ------------
    Psi_mat_t <- as.vector(rowSums(psi_mat_t))
    ## For validated subjects, Psi_kyi = I(Xi=xk, Yi=y) ---------------
    Psi_mat_t[1:n_v] <- 1
    ## -------------------- Update the Psi_kyi for unvalidated subjects
    ## Check that they sum to 1 
    ## data.frame(id = c(seq(1,n_v), rep(seq((n_v+1), n), times = 2*m)), Psi_mat_t) %>% group_by(id) %>% summarize(sumPsi = sum(Psi_mat_t)) %>% select(sumPsi) %>% unique()
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta ---------------------------------------------------
    ### Gradient ------------------------------------------------------
    theta_design_mat <- cbind(int = 1, data_mat[,c(X_val, C)])
    mu <- theta_design_mat %*% prev_theta
    gradient_theta <- matrix(data = c(colSums(Psi_mat_t * c((data_mat[,c(Y_val)]-1 + exp(-mu)/(1+exp(-mu)))) * theta_design_mat)), ncol = 1)
    #if(VERBOSE) print(gradient_theta)
    ### Hessian -------------------------------------------------------
    hessian_theta <- matrix(0, nrow = ncol(theta_design_mat), ncol = ncol(theta_design_mat), byrow = TRUE)
    for (l in 1:ncol(theta_design_mat))
    {
      hessian_theta[l,] <- colSums(c(Psi_mat_t*theta_design_mat[,l])*theta_design_mat*c((exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1)))
    }
    #if(VERBOSE) print(hessian_theta)
    new_theta <- tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta, 
                          error = function(err) {matrix(NA, nrow = nrow(prev_theta))
                          })
    if (TRUE %in% is.na(new_theta))
    {
      warning("Hessian(theta) is not invertible. Iteration stopped before convergence.")
      CONVERGED <- FALSE
      CONVERGED_MSG <- "Hessian theta is not invertible"
      return(list(Coefficients = data.frame(Coefficient = rep(NA, length(prev_theta)), SE = NA),
                  converged = CONVERGED,
                  converged_message = CONVERGED_MSG))
    }
    #if(VERBOSE) print(new_theta)
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta)<TOL
    if(VERBOSE & mean(theta_conv) == 1) print("Theta converged")
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update P(Y*|X*,Y,X,C) ------------------------------------------
    if (NP_YstarDist)
    {
      new_pYstar <- p0
      for (i in 1:nrow(new_pYstar))
      {
        which_Yunval <- as.vector(which(data_mat[-c(1:n_v),Y_unval] == new_pYstar[i,Y_unval]))
        which_Yval <- as.vector(which(data_mat[-c(1:n_v),Y_val] == new_pYstar[i,Y_val]))
        
        # allow for multiple Xs 
        which_Xunval <- as.vector(which(data_mat[-c(1:n_v),X_unval[1]] == new_pYstar[i,X_unval[1]]))
        which_Xval <- as.vector(which(data_mat[-c(1:n_v),X_val[1]] == new_pYstar[i,X_val[1]]))
        if (length(X_val) > 1)
        {
          for (d in 2:length(X_val))
          {
            which_Xunval %<>% intersect(as.vector(which(data_mat[-c(1:n_v),X_unval[d]] == new_pYstar[i,X_unval[d]])))
            which_Xval %<>% as.vector(which(data_mat[-c(1:n_v),X_val[d]] == new_pYstar[i,X_val[d]]))
          }
        }
        new_pYstar[i,"Num_Freq"] <- new_pYstar[i,"Num_Freq"] + sum(Psi_mat_t[intersect(intersect(intersect(which_Xunval, which_Yval), which_Xval), which_Yunval)])
        new_pYstar[i,"Denom_Freq"] <- new_pYstar[i,"Denom_Freq"] + sum(Psi_mat_t[intersect(intersect(which_Xunval, which_Yval), which_Xval)])
      }
      new_pYstar[,"pYstar"] <- new_pYstar[,"Num_Freq"]/new_pYstar[,"Denom_Freq"]
      #if(VERBOSE) print(new_p)
      # Check for convergence -----------------------------------------
      p_conv <- abs(new_pYstar[,"pYstar"] - prev_pYstar[,"pYstar"]) < TOL
      if(VERBOSE & mean(p_conv) == 1) print("{P(Y*|X*,Y,X) coverged")
      new_gamma <- NULL
      gamma_conv <- 1
    } else
    {
      ### Gradient ----------------------------------------------------
      gamma_design_mat <- cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)])
      mu <- gamma_design_mat %*% prev_gamma
      gradient_gamma <- matrix(data = c(colSums(Psi_mat_t * c((data_mat[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
      ### Hessian -----------------------------------------------------
      hessian_gamma <- matrix(0, nrow = ncol(gamma_design_mat), ncol = ncol(gamma_design_mat), byrow = TRUE)
      for (l in 1:ncol(gamma_design_mat))
      {
        hessian_gamma[l,] <- colSums(c(Psi_mat_t*gamma_design_mat[,l])*gamma_design_mat*c((exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1)))
      }
      new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma, 
               error = function(err) {matrix(NA, nrow = nrow(prev_gamma))})
      if (TRUE %in% is.na(new_gamma))
      {
        warning("Hessian(gamma) is not invertible. Iteration stopped before convergence.")
        CONVERGED <- FALSE
        CONVERGED_MSG <- "Hessian gamma is not invertible"
        return(list(Coefficients = data.frame(Coefficient = rep(NA, length(prev_theta)), SE = NA),
                    converged = CONVERGED,
                    converged_message = CONVERGED_MSG))
      }
      p_conv <- 1
      new_pYstar <- NULL
      gamma_conv <- abs(new_gamma - prev_gamma)<TOL
      if(VERBOSE & mean(gamma_conv) == 1) print("Gamma converged")
    }
    ## ------------------------------------------ Update P(Y*|X*,Y,X,C)
    ###################################################################
    ## Update {q_kj} --------------------------------------------------
    ### Update numerators ---------------------------------------------
    new_q_num <- q_val_num + rowsum(psi_mat_t[-c(1:n_v),], group = rep(seq(1,m), each = 2*(n-n_v)))
    ### Update denominators -------------------------------------------
    new_q_denom <- q_val_denom + matrix(data = rep(colSums(psi_mat_t[-c(1:n_v),]), each = m), nrow = m, ncol = sn)
    ### And finally divide them! --------------------------------------
    new_q <- new_q_num/new_q_denom
    #if(VERBOSE) print(new_q)
    ### Check for convergence -----------------------------------------
    q_conv <- abs(new_q - prev_q)<TOL
    #q_conv <- abs(new_q - prev_q)<1E-3
    if(VERBOSE & mean(q_conv) == 1) print("{q_kj} converged")
    ## -------------------------------------------------- Update {q_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(it%%10 == 0 & VERBOSE) print(paste("Iteration", it, "complete."))
    
    all_conv <- c(theta_conv, p_conv, gamma_conv, q_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE

    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_pYstar <- new_pYstar
    prev_gamma <- new_gamma
    prev_q <- new_q
  }
  # ---------------------------------------------- Estimate theta using EM
  if (it >= MAX_ITER & !CONVERGED)
  {
    CONVERGED_MSG <- "MAX_ITER reached without convergence"
    noSE <- TRUE
  }
  
  if (!noSE)
  {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_n <- h_n_scale*n^(-1/2) # perturbation -------------------------------
    ## Calculate pl(theta) -------------------------------------------------
    pl_params <- profile_out(theta = new_theta, 
                             n_v = n_v, n = n,
                             Y_unval=Y_unval, Y_val=Y_val, 
                             X_unval=X_unval, X_val=X_val, 
                             C=C, Bspline=Bspline, 
                             data_mat = data_mat, 
                             p0 = p0, gamma0 = gamma0, q0 = q0, q_val_num = q_val_num, 
                             q_val_denom = q_val_denom, NP_YstarDist = NP_YstarDist)
    
    if(pl_params$converged)
    {
      od_loglik_theta <- observed_data_loglik(n = n, n_v = n_v, 
                                              Y_unval=Y_unval, Y_val=Y_val, 
                                              X_unval=X_unval, X_val=X_val, 
                                              C=C, Bspline=Bspline, 
                                              data_mat = data_mat, 
                                              theta = new_theta, 
                                              p = pl_params$p_at_conv, 
                                              gamma = pl_params$gamma, 
                                              q = pl_params$q_at_conv, 
                                              NP_YstarDist = NP_YstarDist)
    } else
    {
      return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                            SE = NA),
                  converged = CONVERGED,
                  converged_message = "Profile likelihood did not converge. Can't estimate SEs.",
                  NP_YstarDist = NP_YstarDist))
    }
    
    I_theta <- matrix(0, nrow = nrow(new_theta), ncol = nrow(new_theta))
    for (k in 1:ncol(I_theta))
    {
      pert_k <- new_theta; pert_k[k] <- pert_k[k] + h_n
      pl_params <- profile_out(theta = pert_k, 
                               n_v = n_v, n = n,
                               Y_unval=Y_unval, Y_val=Y_val,
                               X_unval=X_unval, X_val=X_val, 
                               C=C, Bspline=Bspline, 
                               data_mat = data_mat, 
                               p0 = p0, gamma0 = gamma0, q0 = q0, q_val_num = q_val_num, q_val_denom = q_val_denom, NP_YstarDist = NP_YstarDist)
      
      if(pl_params$converged)
      {
        od_loglik_pert_k <- observed_data_loglik(n = n, n_v = n_v, 
                                                 Y_unval=Y_unval, Y_val=Y_val, 
                                                 X_unval=X_unval, X_val=X_val, 
                                                 C=C, Bspline=Bspline, 
                                                 data_mat = data_mat, 
                                                 theta = pert_k, 
                                                 p = pl_params$p_at_conv, 
                                                 gamma = pl_params$gamma, 
                                                 q = pl_params$q_at_conv, 
                                                 NP_YstarDist = NP_YstarDist)
      } else
      {
        return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                              SE = NA),
                    converged = CONVERGED,
                    converged_message = "Profile likelihood did not converge. Can't estimate SEs.",
                    NP_YstarDist = NP_YstarDist))
      }
      
      for (l in k:nrow(I_theta))
      {
        pert_l <- new_theta; pert_l[l] <- pert_l[l] + h_n
        pert_both <- pert_l; pert_both[k] <- pert_both[k] + h_n
        
        pl_params <- profile_out(theta = pert_both, 
                                 n_v = n_v, n = n,
                                 Y_unval=Y_unval, Y_val=Y_val, 
                                 X_unval=X_unval, X_val=X_val, 
                                 C=C, Bspline=Bspline, 
                                 data_mat = data_mat, 
                                 p0 = p0, gamma0 = gamma0, q0 = q0, q_val_num = q_val_num, q_val_denom = q_val_denom, NP_YstarDist = NP_YstarDist)
        if (pl_params$converged)
        {
          od_loglik_pert_both <- observed_data_loglik(n_v = n_v, n = n,
                                                      Y_unval=Y_unval, Y_val=Y_val, 
                                                      X_unval=X_unval, X_val=X_val, 
                                                      C=C, Bspline=Bspline, 
                                                      data_mat = data_mat, 
                                                      theta = pert_both, 
                                                      p = pl_params$p_at_conv, 
                                                      gamma = pl_params$gamma, 
                                                      q = pl_params$q_at_conv, 
                                                      NP_YstarDist = NP_YstarDist)
        } else
        {
          return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                                SE = NA),
                      converged = CONVERGED,
                      converged_message = "Profile likelihood did not converge. Can't estimate SEs.",
                      NP_YstarDist = NP_YstarDist))
        }
        
        if (l == k)
        {
          I_theta[l,k] <- od_loglik_pert_both - 2*od_loglik_pert_k + od_loglik_theta
        } else
        {
          pl_params <- profile_out(theta = pert_l, 
                                   n_v = n_v, n = n,
                                   Y_unval=Y_unval, Y_val=Y_val, 
                                   X_unval=X_unval, X_val=X_val, 
                                   C=C, Bspline=Bspline, 
                                   data_mat = data_mat, 
                                   p0 = p0, gamma0 = gamma0, q0 = q0, q_val_num = q_val_num, q_val_denom = q_val_denom, NP_YstarDist = NP_YstarDist)
          
          if (pl_params$converged)
          {
            od_loglik_pert_l <- observed_data_loglik(n_v = n_v, n = n,
                                                     Y_unval=Y_unval, Y_val=Y_val, 
                                                     X_unval=X_unval, X_val=X_val, 
                                                     C=C, Bspline=Bspline, 
                                                     data_mat = data_mat, 
                                                     theta = pert_l, 
                                                     p = pl_params$p_at_conv, 
                                                     gamma = pl_params$gamma, 
                                                     q = pl_params$q_at_conv, 
                                                     NP_YstarDist = NP_YstarDist)
            # symmetry of covariance 
            I_theta[k,l] <- I_theta[l,k] <- od_loglik_pert_both - od_loglik_pert_k - od_loglik_pert_l + od_loglik_theta 
          } else 
          {
            return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                                  SE = NA),
                        converged = CONVERGED,
                        converged_message = "Profile likelihood did not converge. Can't estimate SEs.",
                        NP_YstarDist = NP_YstarDist))
          }
        }
      }
    }
    I_theta <- h_n^(-2) * I_theta
    cov_theta <- tryCatch(expr = -solve(I_theta), 
                          error = function(err) {NA})
    # ------------------------- Estimate Cov(theta) using profile likelihood
    return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                          SE = sqrt(diag(cov_theta))),
                converged = CONVERGED,
                converged_message = CONVERGED_MSG,
                NP_YstarDist = NP_YstarDist))
  } else
  {
    return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                          SE = NA),
                converged = CONVERGED,
                converged_message = CONVERGED_MSG,
                NP_YstarDist = NP_YstarDist))
  }
}
