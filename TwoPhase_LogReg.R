library(magrittr)
library(tidyr)
library(dplyr)

observed_data_loglik <- function(n_v, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data_mat, theta, p, gamma, q, NP_YstarDist)
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
      dplyr::left_join(p) %>% 
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
  log_q <- log(q)
  log_q[log_q == -Inf] <- 0
  log_q %<>% data.frame() %>% magrittr::set_colnames(paste0("log_q_bs", seq(1,sn))) %>% 
    dplyr::mutate(X = unique(data_mat[,X_val]))
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
                       dplyr::left_join(p) %>% 
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
    dplyr::mutate(X = unique(data_mat[,X_val]))
  suppressMessages(unval_data_mat <- data_mat[-c(1:n_v),] %>% data.frame() %>% dplyr::left_join(q) %>% data.matrix())
  ## --------------------------------------------------------------- Merge in q_kj
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x q_kj ---------------------------
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(c(pY_X*phat)*unval_data_mat[,Bspline] * unval_data_mat[,paste0("q_bs", seq(1,sn))])
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}

profile_out <- function(theta, n_v, n, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data_mat, p0, gamma0, q0, q_val_num, q_val_denom, NP_YstarDist, VERBOSE = FALSE, TOL = 1E-4, MAX_ITER = 200)
{
  sn <- ncol(q0)
  m <- nrow(q0)
  
  prev_pYstar <- p0
  prev_gamma <- gamma0
  prev_q <- q0
  
  CONVERGED <- FALSE
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
                         dplyr::left_join(data.frame(prev_pYstar, row.names = NULL)) %>% 
                         dplyr::pull(pYstar))
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
    prev_logq_stacked <- log(prev_q_stacked)
    prev_logq_stacked[prev_logq_stacked == -Inf] <- 0
    ### ---------------------------------------------------------- q_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)B(X*)log(q) ----------------------------
    psi_mat_t_num <- pY_X * phat * data_mat[,Bspline] * prev_logq_stacked
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
        which_Xunval <- as.vector(which(data_mat[-c(1:n_v),X_unval] == new_pYstar[i,X_unval]))
        which_Yval <- as.vector(which(data_mat[-c(1:n_v),Y_val] == new_pYstar[i,Y_val]))
        which_Xval <- as.vector(which(data_mat[-c(1:n_v),X_val] == new_pYstar[i,X_val]))
        
        new_pYstar$Num_Freq[i] <- new_pYstar$Num_Freq[i] + sum(Psi_mat_t[intersect(intersect(intersect(which_Xunval, which_Yval), which_Xval), which_Yunval)])
        new_pYstar$Denom_Freq[i] <- new_pYstar$Denom_Freq[i] + sum(Psi_mat_t[intersect(intersect(which_Xunval, which_Yval), which_Xval)])
      }
      new_pYstar %<>% dplyr::mutate(pYstar = Num_Freq/Denom_Freq)
      #if(VERBOSE) print(new_p)
      # Check for convergence -----------------------------------------
      p_conv <- abs(new_pYstar$pYstar - prev_pYstar$pYstar)<TOL
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
      if (abs(det(hessian_gamma)) < 1e-08)
      {
        new_gamma <- matrix(NA, nrow = nrow(prev_gamma))
        HESSIAN_INVERTIBLE <- FALSE
        break
      } else 
      {
        new_gamma <- prev_gamma - solve(hessian_gamma) %*% gradient_gamma
        HESSIAN_INVERTIBLE <- TRUE
      }
      #if(VERBOSE) print(new_gamma)
      # Check for convergence -----------------------------------------
      p_conv <- 1
      new_pYstar <- NULL
      gamma_conv <- abs(new_gamma - prev_gamma)<TOL
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
  return(list("psi_mat_at_conv" = psi_mat_t,
              "p_at_conv" = new_pYstar, 
              "gamma_at_conv" = new_gamma,
              "q_at_conv" = new_q))
}

TwoPhase_LogReg <- function(Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data, h_n_scale = 1, VERBOSE = FALSE, TOL = 1E-4, MAX_ITER = 200)
{
  n <- nrow(data)
  V <- !is.na(data[,X_val]) & !is.na(data[Y_val]) # validation indicators
  n_v <- sum(V)
  id <- seq(1,n)
  
  # Reorder so that the n_v validated subjects are first ------------
  data <- data[order(as.numeric(V), decreasing = TRUE),]
  
  # If unvalidated variable was left blank, assume error-free -------
  if (is.null(Y_unval)) Y_unval <- Y_val
  if (is.null(X_unval)) X_unval <- X_val
  
  Ystar <- data[,Y_unval]
  Y <- data[,Y_val]
  Xstar <- data[,X_unval]
  X <- data[,X_val]
  #if(!is.null(C)) Covar <- data[,C] else Covar <- rep(0,nrow(data))
  
  # Add the B spline basis ------------------------------------------
  B <- data[,Bspline]
  sn <- ncol(B)
  # ------------------------------------------ Add the B spline basis
  
  # Add m for distinct values of X ----------------------------------
  x_k <- unique(data.frame(X)[1:n_v,])
  m <- nrow(data.frame(x_k))
  # ---------------------------------- Add m for distinct values of X 
  
  p <- length(X_val) + length(C) # number of coefficients in the analysis model (beta_1, ..., beta_p)
  prev_theta <- matrix(0, nrow = (p+1), ncol = 1)
  
  if (0 %in% data.frame(table(data[c(1:n_v), c(X_unval, Y_val, X_val, C)]))$Freq)
  {
    warning("Logistic regression will be used to estimate P(Y*|X*,Y,X,C).")
    NP_YstarDist <- FALSE
    prev_gamma <- gamma0 <- matrix(0, nrow = (length(Y_val) + length(X_unval) + length(X_val) + length(C) + 1), ncol = 1)
    p0 <- prev_pYstar <- NULL
  } else
  {
    NP_YstarDist <- TRUE
    tab_all <- data.frame(table(data[c(1:n_v), c(Y_unval, X_unval, Y_val, X_val, C)]))
    colnames(tab_all)[length(colnames(tab_all))] <- "Num_Freq"
    tab_exclYstar <- data.frame(table(data[c(1:n_v), c(X_unval, Y_val, X_val, C)]))
    colnames(tab_exclYstar)[length(colnames(tab_exclYstar))] <- "Denom_Freq"
    suppressMessages(p0 <- prev_pYstar <- dplyr::inner_join(tab_all, tab_exclYstar) %>% 
      dplyr::mutate(pYstar = Num_Freq/Denom_Freq, 
                    Ystar = as.numeric(as.character(Ystar)), Y = as.numeric(as.character(Y)),
                    Xstar = as.numeric(as.character(Xstar)), X = as.numeric(as.character(X))))
    gamma0 <- prev_gamma <- NULL
  }
  
  # Replace initial q w/ Bspline version  ------------ 
  prev_q <- q_val_num <- q_val_denom <- matrix(0, nrow = m, ncol = sn)
  for (j in 1:sn)
  {
    for (k in 1:m)
    {
      prev_q[k,j] <- q_val_num[k,j] <- sum(as.numeric(X[1:n_v] == x_k[k])*B[V,j])
    }
    q_val_denom[,j] <- sum(prev_q[,j]) 
    prev_q[,j] <- prev_q[,j]/sum(prev_q[,j])
  }
  colnames(prev_q) <- paste0("q_bs", seq(1,sn))
  q0 <- prev_q
  # ------------ Replace initial q w/ Bspline version 
  
  # Add the Bspline to the psi_mat_t ----------------
  ## For validated subjects, keep only Xi/Yi --------
  ## (So one row per person) ------------------------
  ## For unvalidated subjects, create a row for -----
  ## each observed (xk,y) ---------------------------
  if(!is.null(C))
  {
    data_mat <- data.matrix(cbind(data[c(1:n_v),c(X_val, Y_val, X_unval, Y_unval, C)], B[c(1:n_v),]))
    suppressWarnings(data_mat %<>% rbind(data.matrix(cbind(X = rep(x_k, each = 2*(n-n_v)),
                                                           Y = rep(rep(c(0,1), each = (n-n_v)), times = m),
                                                           data[-c(1:n_v),c(X_unval, Y_unval, C)], 
                                                           B[-c(1:n_v),]))))
  } else 
  {
      data_mat <- data.matrix(cbind(X, Y, Xstar, Ystar, B)[1:n_v,])
      suppressWarnings(data_mat %<>% rbind(data.matrix(cbind(X = rep(x_k, each = 2*(n-n_v)),
                                                              Y = rep(rep(c(0,1), each = (n-n_v)), times = m),
                                                              data[-c(1:n_v),c(X_unval, Y_unval)], 
                                                              B[-c(1:n_v),])
      )))
  }
  # ---------------- Add the Bspline to the psi_mat_t
  
  CONVERGED <- FALSE
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
      suppressMessages(phat <- data_mat %>% data.frame(row.names = NULL) %>% dplyr::left_join(data.frame(prev_pYstar, row.names = NULL)) %>% dplyr::pull(pYstar))
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
    prev_logq_stacked <- log(prev_q_stacked)
    prev_logq_stacked[prev_logq_stacked == -Inf] <- 0
    ### ---------------------------------------------------------- q_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)B(X*)log(q) ----------------------------
    psi_mat_t_num <- pY_X * phat * data_mat[,Bspline] * prev_logq_stacked
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
    ## data.frame(id = c(seq(1,n_v), rep(seq((n_v+1), n), times = 2*m)), Psi_mat_t) %>% group_by(id) %>% summarize(sumPsi = sum(Psi_mat_t)) %>% arrange(desc(sumPsi))
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
    if (abs(det(hessian_theta)) < 1e-08)
    {
      new_theta <- matrix(NA, nrow = nrow(prev_theta))
      HESSIAN_INVERTIBLE <- FALSE
      break
    } else 
    {
      new_theta <- prev_theta - solve(hessian_theta) %*% gradient_theta
      HESSIAN_INVERTIBLE <- TRUE
    }
    if(VERBOSE) print(new_theta)
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta)<TOL
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update P(Y*|X*,Y,X,C) ------------------------------------------
    if (NP_YstarDist)
    {
      new_pYstar <- p0
      for (i in 1:nrow(new_pYstar))
      {
        which_Yunval <- as.vector(which(data_mat[-c(1:n_v),Y_unval] == new_pYstar[i,Y_unval]))
        which_Xunval <- as.vector(which(data_mat[-c(1:n_v),X_unval] == new_pYstar[i,X_unval]))
        which_Yval <- as.vector(which(data_mat[-c(1:n_v),Y_val] == new_pYstar[i,Y_val]))
        which_Xval <- as.vector(which(data_mat[-c(1:n_v),X_val] == new_pYstar[i,X_val]))
        
        new_pYstar$Num_Freq[i] <- new_pYstar$Num_Freq[i] + sum(Psi_mat_t[intersect(intersect(intersect(which_Xunval, which_Yval), which_Xval), which_Yunval)])
        new_pYstar$Denom_Freq[i] <- new_pYstar$Denom_Freq[i] + sum(Psi_mat_t[intersect(intersect(which_Xunval, which_Yval), which_Xval)])
      }
      new_pYstar %<>% dplyr::mutate(pYstar = Num_Freq/Denom_Freq)
      #if(VERBOSE) print(new_p)
      # Check for convergence -----------------------------------------
      p_conv <- abs(new_pYstar$pYstar - prev_pYstar$pYstar)<TOL
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
      if (abs(det(hessian_gamma)) < 1e-08)
      {
        new_gamma <- matrix(NA, nrow = nrow(prev_gamma))
        HESSIAN_INVERTIBLE <- FALSE
        break
      } else 
      {
        new_gamma <- prev_gamma - solve(hessian_gamma) %*% gradient_gamma
        HESSIAN_INVERTIBLE <- TRUE
      }
      #if(VERBOSE) print(new_gamma)
      # Check for convergence -----------------------------------------
      p_conv <- 1
      new_pYstar <- NULL
      gamma_conv <- abs(new_gamma - prev_gamma)<TOL
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
    ## -------------------------------------------------- Update {q_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(VERBOSE) print(paste("Iteration", it, "complete."))
    
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
  # Estimate Cov(theta) using profile likelihood -------------------------
  h_n <- h_n_scale*n^(-1/2) # perturbation -------------------------------
  ## Calculate pl(theta) -------------------------------------------------
  pl_params <- profile_out(theta = new_theta, 
                           n_v = n_v, n = n,
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           data_mat = data_mat, 
                           p0 = p0, gamma0 = gamma0, q0 = q0, q_val_num = q_val_num, q_val_denom = q_val_denom, NP_YstarDist = NP_YstarDist)
  
  od_loglik_theta <- observed_data_loglik(n_v = n_v, 
                                          Y_unval=Y_unval, Y_val=Y_val, 
                                          X_unval=X_unval, X_val=X_val, 
                                          C=C, Bspline=Bspline, 
                                          data_mat = data_mat, 
                                          theta = new_theta, 
                                          p = pl_params$p_at_conv, 
                                          gamma = pl_params$gamma, 
                                          q = pl_params$q_at_conv, 
                                          NP_YstarDist = NP_YstarDist)
  
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
    od_loglik_pert_k <- observed_data_loglik(n_v = n_v, 
                                             Y_unval=Y_unval, Y_val=Y_val, 
                                             X_unval=X_unval, X_val=X_val, 
                                             C=C, Bspline=Bspline, 
                                             data_mat = data_mat, 
                                             theta = pert_k, 
                                             p = pl_params$p_at_conv, 
                                             gamma = pl_params$gamma, 
                                             q = pl_params$q_at_conv, 
                                             NP_YstarDist = NP_YstarDist)
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
      od_loglik_pert_both <- observed_data_loglik(n_v = n_v, 
                                                  Y_unval=Y_unval, Y_val=Y_val, 
                                                  X_unval=X_unval, X_val=X_val, 
                                                  C=C, Bspline=Bspline, 
                                                  data_mat = data_mat, 
                                                  theta = pert_both, 
                                                  p = pl_params$p_at_conv, 
                                                  gamma = pl_params$gamma, 
                                                  q = pl_params$q_at_conv, 
                                                  NP_YstarDist = NP_YstarDist)
      
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
        od_loglik_pert_l <- observed_data_loglik(n_v = n_v, 
                                                 Y_unval=Y_unval, Y_val=Y_val, 
                                                 X_unval=X_unval, X_val=X_val, 
                                                 C=C, Bspline=Bspline, 
                                                 data_mat = data_mat, 
                                                 theta = pert_l, 
                                                 p = pl_params$p_at_conv, 
                                                 gamma = pl_params$gamma, 
                                                 q = pl_params$q_at_conv, 
                                                 NP_YstarDist = NP_YstarDist)
        I_theta[k,l] <- I_theta[l,k] <- od_loglik_pert_both - od_loglik_pert_k - od_loglik_pert_l + od_loglik_theta # symmetry of covariance 
      }
    }
  }
  I_theta <- h_n^(-2) * I_theta
  cov_theta <- -solve(I_theta)
  # ------------------------- Estimate Cov(theta) using profile likelihood
  
  return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                        SE = sqrt(diag(cov_theta))),
              Est_YstarDist_Nonparametrically = NP_YstarDist, 
              converged = CONVERGED))
}

