library(magrittr); library(tidyr); library(dplyr)

profile_out <- function(theta, n_v, n, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data_mat, gamma0, p0, p_val_num, p_val_denom, VERBOSE = FALSE, TOL = 1E-4, MAX_ITER = 1000, useGLM = FALSE)
{
  sn <- ncol(p0)
  m <- nrow(p0)
  
  prev_theta <- theta
  prev_gamma <- gamma0
  prev_p <- p0
  
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  
  CONVERGED <- FALSE
  it <- 1
  
  # Estimate gamma/p using EM ----------------------------------------------
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
    pYstar <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)]) %*% prev_gamma)))
    pYstar <- ifelse(as.vector(data_mat[,Y_unval]) == 0, 1-pYstar, pYstar)
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    suppressMessages(pX <- data.frame(data_mat) %>% dplyr::left_join(data.frame(k = 1:m, prev_p)) %>% dplyr::arrange(orig_row))
    pX <- pX[,paste0("p_k", seq(1,sn))]
    # make it 0 for validated to avoid confusion
    pX[c(1:n_v),] <- 0
    ### ---------------------------------------------------------- p_kj
    ###################################################################
    ### psi_kjyi ------------------------------------------------------
    psi_t_num <- pY_X * pYstar * data_mat[,Bspline] * pX
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = data_mat[,"id"], reorder = TRUE)
    psi_t_num_sumover_k_y_j <- rowSums(psi_t_num_sumover_k_y)
    suppressMessages(psi_t_denom_stacked <- data.frame(data_mat) %>% dplyr::left_join(data.frame(id = 1:n, psi_denom = psi_t_num_sumover_k_y_j)) %>% dplyr::select(psi_denom) %>% data.matrix())
    psi_t_denom_stacked[c(1:n_v)] <- 1
    psi_t_denom_wide <- matrix(rep(psi_t_denom_stacked, times = sn), nrow = nrow(psi_t_denom_stacked), ncol = sn, byrow = FALSE)
    psi_t <- psi_t_num/psi_t_denom_wide
    ### ------------------------------------------------------ psi_kjyi
    ###################################################################
    ## Update the w_kyi for unvalidated subjects ----------------------
    w_t <- rowSums(psi_t)
    w_t[1:n_v] <- 1
    ki <- paste0("k", data_mat[,"k"], "_id", data_mat[,"id"])
    u_t <- rowsum(psi_t, group = ki, reorder = FALSE)
    ## ---------------------- Update the u_kji for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update gamma ---------------------------------------------------
    if(!useGLM)
    {
      ### Gradient ------------------------------------------------------
      gamma_design_mat <- cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)])
      mu <- gamma_design_mat %*% prev_gamma
      gradient_gamma <- matrix(data = c(colSums(w_t * c((data_mat[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
      ### ------------------------------------------------------ Gradient
      ### Hessian -------------------------------------------------------
      hessian_gamma <- matrix(0, nrow = ncol(gamma_design_mat), ncol = ncol(gamma_design_mat), byrow = TRUE)
      post_multiply <- c(w_t*(exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))
      for (l in 1:ncol(gamma_design_mat))
      {
        hessian_gamma[l,] <- colSums(c(gamma_design_mat[,l])*gamma_design_mat*post_multiply)
      }
      
      new_gamma <- tryCatch(expr = prev_gamma - newton_step_scale *solve(hessian_gamma) %*% gradient_gamma,
                            error = function(err) {matrix(NA, nrow = nrow(prev_gamma))
                            })
      if (TRUE %in% is.na(new_gamma))
      {
        # CONVERGED <- FALSE
        # CONVERGED_MSG <- "Hessian gamma is nonsingular"
        # return(list(Coefficients = data.frame(Coefficient = rep(NA, length(prev_theta)), 
        #                                       SE = NA),
        #             converged = CONVERGED,
        #             converged_msg = CONVERGED_MSG))
        useGLM <- TRUE
      }
    }
    
    if(useGLM)
    {
      suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(data_mat), weights = w_t)$coefficients, ncol = 1))
    }
    
    # ## Update gamma ---------------------------------------------------
    # if(useGLM)
    # {
    #   suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(data_mat), weights = w_t)$coefficients, ncol = 1))
    # }
    # else
    # {
    #   ### Gradient ------------------------------------------------------
    #   gamma_design_mat <- cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)])
    #   mu <- gamma_design_mat %*% prev_gamma
    #   gradient_gamma <- matrix(data = c(colSums(w_t * c((data_mat[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
    #   ### ------------------------------------------------------ Gradient
    #   ### Hessian -------------------------------------------------------
    #   hessian_gamma <- matrix(0, nrow = ncol(gamma_design_mat), ncol = ncol(gamma_design_mat), byrow = TRUE)
    #   post_multiply <- c(w_t*(exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))
    #   for (l in 1:ncol(gamma_design_mat))
    #   {
    #     hessian_gamma[l,] <- colSums(c(gamma_design_mat[,l])*gamma_design_mat*post_multiply)
    #   }
    #   
    #   new_gamma <- tryCatch(expr = prev_gamma - newton_step_scale*solve(hessian_gamma) %*% gradient_gamma,
    #                         error = function(err) {matrix(NA, nrow = nrow(prev_gamma))
    #                         })
    # }
    # 
    # if (TRUE %in% is.na(new_gamma))
    # {
    #   CONVERGED <- FALSE
    #   CONVERGED_MSG <- "Hessian gamma is nonsingular"
    #   return(list(Coefficients = data.frame(Coefficient = rep(NA, length(prev_gamma)), SE = NA),
    #               converged = CONVERGED,
    #               converged_msg = CONVERGED_MSG))
    # }
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma)<TOL
    ## ------------------------------------------ Update gamma
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    new_p_num <- p_val_num + rowsum(u_t[-c(1:n_v),], group = rep(seq(1,m), times = (n-n_v)), reorder = TRUE)
    ### Update denominators by summing new_p_num w/i each col ---------
    new_p_denom <- matrix(rep(colSums(new_p_num), times = m), nrow = m, ncol = sn, byrow = TRUE)
    new_p <- new_p_num/new_p_denom
    colnames(new_p) <- colnames(prev_p)
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p)<TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(VERBOSE & it%%25 == 0) print(paste("Iteration", it, "complete."))
    
    all_conv <- c(gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_gamma <- new_gamma
    prev_p <- new_p 
  }
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "gamma_at_conv" = new_gamma,
              "p_at_conv" = new_p))
}

observed_data_loglik <- function(n_v, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, data_mat, theta, gamma, p)
{
  sn <- ncol(p)
  m <- nrow(p)
  
  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[c(1:n_v),c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(data_mat[c(1:n_v),c(Y_val)]) == 0, 1-pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over log[P(Yi*|Xi*,Yi,Xi)] -----------------------------------------------
  pYstar <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[c(1:n_v),c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar <- ifelse(as.vector(data_mat[c(1:n_v),Y_unval]) == 0, 1-pYstar, pYstar)
  return_loglik <- return_loglik + sum(log(pYstar))
  ## ----------------------------------------------- Sum over log[P(Yi*|Xi*,Yi,Xi)]
  #################################################################################
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
  sn <- ncol(p)
  suppressMessages(pX <- data.frame(data_mat[c(1:n_v),]) %>% dplyr::left_join(data.frame(k = 1:m, p)) %>% dplyr::select(paste0("p_k", seq(1,sn))))
  log_pX <- log(pX)
  log_pX[log_pX == -Inf] <- 0
  return_loglik <- return_loglik + sum(data_mat[c(1:n_v),Bspline] * log_pX)
  ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  #################################################################################
  # -------------------------------------------------------- For validated subjects
  
  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[-c(1:n_v),c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(data_mat[-c(1:n_v),c(Y_val)]) == 0, 1-pY_X, pY_X)
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate P(Yi*|Xi*,y,xk) for all (y,xk) ------------------------------------
  pYstar <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[-c(1:n_v),c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar <- ifelse(as.vector(data_mat[-c(1:n_v),Y_unval]) == 0, 1-pYstar, pYstar)
  ## ------------------------------------ Calculate P(Yi*|Xi*,y,xk) for all (y,xk)
  ################################################################################
  ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
  suppressMessages(pX <- data.frame(data_mat[-c(1:n_v),]) %>% 
                     dplyr::left_join(data.frame(k = 1:m, p)))
  ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x q_kj ---------------------------
  person_sum <- rowsum(c(pY_X*pYstar)*pX[,Bspline]*pX[,paste0("p_k", seq(1,sn))], group = data_mat[-c(1:n_v), "id"])
  log_person_sum <- log(person_sum)
  log_person_sum[log_person_sum == -Inf] <- 0
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}

TwoPhase_LogReg <- function(Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, ValidationIndicator = NULL, Bspline=NULL, data, h_n_scale = 1, newton_step_scale = 1, noSE=FALSE, useGLM = FALSE, VERBOSE = FALSE, TOL = 1E-6, MAX_ITER = 1000, save_all_iter = FALSE)
{
  n <- nrow(data)
  #V <- data[,ValidationIndicator]
  n_v <- sum(data[,ValidationIndicator])
  
  # Reorder so that the n_v validated subjects are first ------------
  data <- data[order(as.numeric(data[,ValidationIndicator]), decreasing = TRUE),]
  
  # If unvalidated variable was left blank, assume error-free -------
  ## Need to write simplification here 

  # Add the B spline basis ------------------------------------------
  B <- data[,Bspline]
  sn <- ncol(B)
  if(0 %in% colSums(B[c(1:n_v),]))
  {
    return(list(Coefficients = data.frame(Coefficient = NA, 
                                          SE = NA),
                converged = FALSE,
                converged_msg = "Empty sieve in validated data."))
  }
  # ------------------------------------------ Add the B spline basis
  
  # Add m for distinct values of X ----------------------------------
  x_k <- data.frame(unique(data[1:n_v,c(X_val)]))
  x_k <- data.frame(x_k[order(x_k[,1]),])
  colnames(x_k) <- c(X_val)
  m <- nrow(x_k)
  suppressMessages(data %<>% dplyr::left_join(data.frame(x_k, k = 1:m)))
  # ---------------------------------- Add m for distinct values of X 
  
  num_pred <- length(X_val) + length(C) # number of coefficients in the analysis model (beta_1, ..., beta_p)
  prev_theta <- theta0 <- matrix(0, nrow = (num_pred+1), ncol = 1)
  theta_formula <- as.formula(paste0(Y_val, "~", paste(c(X_val, C), collapse = "+")))
  
  prev_gamma <- gamma0 <- matrix(0, nrow = (length(Y_val) + length(X_unval) + num_pred + 1), ncol = 1)
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  
  # Initialize B-spline coefficients {p_kj}  ------------ 
  ## Numerators sum B(Xi*) over k = 1,...,m -------------
  p_val_num <- rowsum(data[c(1:n_v), Bspline], group = data[c(1:n_v),"k"], reorder = TRUE)
  p_val_denom <- matrix(rep(colSums(p_val_num), times = m), nrow = m, ncol = sn, byrow = TRUE)
  p0 <- p_val_num/p_val_denom
  colnames(p0) <- paste0("p_k", seq(1,sn))
  rownames(p0) <- paste0("p_", seq(1,m), "j")
  prev_p <- p0
  # Since reorder = TRUE, the rows of prev_p will be -----
  # k = 1, ..., m ----------------------------------------
  # -------------- Initialize B-spline coefficients {p_kj}
  
  ## For validated subjects, keep only Xi/Yi --------
  ## (one row per person) ---------------------------
  data_mat <- data.matrix(cbind(data[c(1:n_v), c(X_val, Y_val, X_unval, Y_unval, C, Bspline, "k")],id = 1:n_v))
  
  ## For unvalidated subjects, save Xi*, Yi*, Ci, B(Xi*)
  data_mat_unval <- data[-c(1:n_v), c(X_unval, Y_unval, C, Bspline)] %>% 
    dplyr::mutate(unval = TRUE, id = (n_v+1):n)
  
  ## Could have observed Y = 0 or Y = 1
  suppressMessages(data_mat_unval %<>% dplyr::full_join(data.frame(unval = TRUE, unobs_Y = c(0,1))))
  colnames(data_mat_unval)[length(colnames(data_mat_unval))] <- Y_val
  
  ## Could have observed X=x1, ..., X=xk
  suppressMessages(data_mat_unval %<>% dplyr::full_join(data.frame(unval = TRUE, x_k, k = 1:m)))
  
  data_mat_unval <- data_mat_unval[,c(X_val, Y_val, X_unval, Y_unval, C, Bspline, "k", "id")]
  
  ## Stack validated on top of unvalidated 
  data_mat %<>% rbind(data_mat_unval) %>% dplyr::mutate(orig_row = 1:n()) %>% data.matrix()

  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  useGLM_theta <- useGLM_gamma <- useGLM
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
    pYstar <- 1/(1 + exp(-as.numeric(cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)]) %*% prev_gamma)))
    pYstar <- ifelse(as.vector(data_mat[,Y_unval]) == 0, 1-pYstar, pYstar)
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    suppressMessages(pX <- data.frame(data_mat) %>% dplyr::left_join(data.frame(k = 1:m, prev_p)) %>% dplyr::arrange(orig_row))
    pX <- pX[,paste0("p_k", seq(1,sn))]
    # make it 0 for validated to avoid confusion
    pX[c(1:n_v),] <- 0
    ### ---------------------------------------------------------- p_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)p_kjB(X*) ------------------------------
    psi_t_num <- pY_X * pYstar * data_mat[,Bspline] * pX
    ### Update denominator --------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ------------------
    #### reorder = TRUE returns them in ascending order of i ----------
    #### (rather than in order of encounter) --------------------------
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = data_mat[,"id"], reorder = TRUE)
    #### Then sum over the sn splines ---------------------------------
    #### Same ordering as psi_t_num_sumover_k_y, just only 1 column ---
    psi_t_num_sumover_k_y_j <- rowSums(psi_t_num_sumover_k_y)
    #### Since psi_t_num is (nv + (n-nv)*2*m) x sn, same as data_mat --
    #### And in the same order, use a left_join() to assign person- ---
    #### specific denominators 
    suppressMessages(psi_t_denom_stacked <- data.frame(data_mat) %>% dplyr::left_join(data.frame(id = 1:n, psi_denom = psi_t_num_sumover_k_y_j)) %>% dplyr::select(psi_denom) %>% data.matrix())
    #### For validated, make sumPsi_over_xk_y_j = 1 rather than 0 -----
    #### to avoid math issues with dividing by 0
    psi_t_denom_stacked[c(1:n_v)] <- 1
    #psi_t_denom_stacked <- ifelse(psi_t_denom_stacked == 0, 1, psi_t_denom_stacked)
    #### Cbind psi_denom_stacked together sn times so that each col ---
    #### of psi_t_num can be divided by it ----------------------------
    psi_t_denom_wide <- matrix(rep(psi_t_denom_stacked, times = sn), nrow = nrow(psi_t_denom_stacked), ncol = sn, byrow = FALSE)
    ### And divide them! ----------------------------------------------
    psi_t <- psi_t_num/psi_t_denom_wide
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
    ## Update the w_kyi for unvalidated subjects ----------------------
    ## by summing across the splines/ columns of psi_t ----------------
    ## w_t is ordered by i = 1, ..., n --------------------------------
    ## and w/i Vi = 0 by y = 0/1 and then by k = 1, ..., m ------------
    ## (same ordering as in data_mat)
    w_t <- rowSums(psi_t)
    ## ---------------------- Update the w_kyi for unvalidated subjects
    ## For validated subjects, w_kyi = I(Xi=xk, Yi=y) ---------------
    w_t[1:n_v] <- 1
    ## Update the u_kji for unvalidated subjects ----------------------
    ## by summing over Y = 0/1 w/i each i, k --------------------------
    ## do this by concatenating k & id columns ------------------------
    ki <- paste0("k", data_mat[,"k"], "_id", data_mat[,"id"])
    ## Set reorder = FALSE so that u_t will be returned in the order --
    ## in which the ki were encountered -------------------------------
    ## (i.e. i = 1,..., n then k = 1, ..., m) -------------------------
    u_t <- rowsum(psi_t, group = ki, reorder = FALSE)
    ## ---------------------- Update the u_kji for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta ---------------------------------------------------
    if(!useGLM_theta)
    {
      ### Gradient ------------------------------------------------------
      theta_design_mat <- cbind(int = 1, data_mat[,c(X_val, C)])
      mu <- theta_design_mat %*% prev_theta
      gradient_theta <- matrix(data = c(colSums(w_t * c((data_mat[,c(Y_val)]-1 + exp(-mu)/(1+exp(-mu)))) * theta_design_mat)), ncol = 1)
      ### ------------------------------------------------------ Gradient
      ### Hessian -------------------------------------------------------
      hessian_theta <- matrix(0, nrow = ncol(theta_design_mat), ncol = ncol(theta_design_mat), byrow = TRUE)
      post_multiply <- c((exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))
      for (l in 1:ncol(theta_design_mat))
      {
        hessian_theta[l,] <- colSums(c(w_t*theta_design_mat[,l])*theta_design_mat*post_multiply)
      }
      
      new_theta <- tryCatch(expr = prev_theta - newton_step_scale*solve(hessian_theta) %*% gradient_theta,
                            error = function(err) {matrix(NA, nrow = nrow(prev_theta))
                            })
      
      if (TRUE %in% is.na(new_theta))
      {
        # CONVERGED <- FALSE
        # CONVERGED_MSG <- "Hessian theta is nonsingular."
        # return(list(Coefficients = data.frame(Coefficient = rep(NA, length(prev_theta)), SE = NA),
        #             converged = CONVERGED,
        #             converged_msg = CONVERGED_MSG))
        useGLM_theta <- TRUE
      }
    }
    
    if(useGLM_theta)
    {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(data_mat), weights = w_t)$coefficients, ncol = 1))
    }
    if(VERBOSE) print(new_theta)
    ### ------------------------------------------------------- Hessian
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta)<TOL
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update gamma ---------------------------------------------------
    if(!useGLM_gamma)
    {
      ### Gradient ------------------------------------------------------
      gamma_design_mat <- cbind(int = 1, data_mat[,c(X_unval, Y_val, X_val, C)])
      mu <- gamma_design_mat %*% prev_gamma
      gradient_gamma <- matrix(data = c(colSums(w_t * c((data_mat[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
      ### ------------------------------------------------------ Gradient
      ### Hessian -------------------------------------------------------
      hessian_gamma <- matrix(0, nrow = ncol(gamma_design_mat), ncol = ncol(gamma_design_mat), byrow = TRUE)
      post_multiply <- c(w_t*(exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))
      for (l in 1:ncol(gamma_design_mat))
      {
        hessian_gamma[l,] <- colSums(c(gamma_design_mat[,l])*gamma_design_mat*post_multiply)
      }
      
      new_gamma <- tryCatch(expr = prev_gamma - newton_step_scale *solve(hessian_gamma) %*% gradient_gamma,
                            error = function(err) {matrix(NA, nrow = nrow(prev_gamma))
                            })
      if (TRUE %in% is.na(new_gamma))
      {
        # CONVERGED <- FALSE
        # CONVERGED_MSG <- "Hessian gamma is nonsingular"
        # return(list(Coefficients = data.frame(Coefficient = rep(NA, length(prev_theta)), 
        #                                       SE = NA),
        #             converged = CONVERGED,
        #             converged_msg = CONVERGED_MSG))
        useGLM_gamma <- TRUE
      }
    }
    
    if(useGLM_gamma)
    {
      suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(data_mat), weights = w_t)$coefficients, ncol = 1))
    }
    if(VERBOSE) print(new_gamma)
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma)<TOL/0.01
    ## ------------------------------------------ Update gamma
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., n -----------
    ### Since reorder = TRUE, rows will be in order k = 1, ..., m -----
    ### Recall: u_t is ordered by i = 1,..., n then k = 1, ..., m -----
    ### So, group = rep(seq(1,m), times = (n-n_v)) (excl. Vi = 1) -----
    
    new_p_num <- p_val_num + rowsum(u_t[-c(1:n_v),], group = rep(seq(1,m), times = (n-n_v)), reorder = TRUE)
    ### Update denominators by summing new_p_num w/i each col ---------
    new_p_denom <- matrix(rep(colSums(new_p_num), times = m), nrow = m, ncol = sn, byrow = TRUE)
    new_p <- new_p_num/new_p_denom
    colnames(new_p) <- colnames(prev_p)
    if(VERBOSE) print(new_p[1,])
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p)<TOL/0.01
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(VERBOSE & it%%25 == 0) print(paste("Iteration", it, "complete."))
    
    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p 
  }
  
  if(!CONVERGED & it > MAX_ITER) 
  {
    CONVERGED_MSG = "MAX_ITER reached"
    new_theta <- matrix(NA, nrow = nrow(prev_theta))
  }
  if(CONVERGED) CONVERGED_MSG <- "Converged"
  
  # ---------------------------------------------- Estimate theta using EM
  if(noSE | !CONVERGED)
  {
      return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                            SE = NA),
                  converged = CONVERGED,
                  converged_msg = CONVERGED_MSG))
  } else
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
                             gamma0 = gamma0, p0 = p0, p_val_num = p_val_num, p_val_denom = p_val_denom, useGLM = useGLM)
    od_loglik_theta <- observed_data_loglik(n_v = n_v, 
                                            Y_unval=Y_unval, Y_val=Y_val, 
                                            X_unval=X_unval, X_val=X_val, 
                                            C=C, Bspline=Bspline, 
                                            data_mat = data_mat, 
                                            theta = new_theta,
                                            gamma = pl_params$gamma, 
                                            p = pl_params$p_at_conv)
    
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
                               gamma0 = gamma0, 
                               p0 = p0, p_val_num = p_val_num, p_val_denom = p_val_denom, useGLM = useGLM)
      od_loglik_pert_k <- observed_data_loglik(n_v = n_v, 
                                               Y_unval=Y_unval, Y_val=Y_val, 
                                               X_unval=X_unval, X_val=X_val, 
                                               C=C, Bspline=Bspline, 
                                               data_mat = data_mat, 
                                               theta = pert_k,
                                               gamma = pl_params$gamma, 
                                               p = pl_params$p_at_conv)
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
                                 gamma0 = gamma0, 
                                 p0 = p0, p_val_num = p_val_num, p_val_denom = p_val_denom, useGLM = useGLM)
        od_loglik_pert_both <- observed_data_loglik(n_v = n_v, 
                                                    Y_unval=Y_unval, Y_val=Y_val, 
                                                    X_unval=X_unval, X_val=X_val, 
                                                    C=C, Bspline=Bspline, 
                                                    data_mat = data_mat, 
                                                    theta = pert_both, 
                                                    gamma = pl_params$gamma, 
                                                    p = pl_params$p_at_conv)
        
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
                                   gamma0 = gamma0, 
                                   p0 = p0, p_val_num = p_val_num, p_val_denom = p_val_denom, useGLM = useGLM)
          od_loglik_pert_l <- observed_data_loglik(n_v = n_v, 
                                                   Y_unval=Y_unval, Y_val=Y_val, 
                                                   X_unval=X_unval, X_val=X_val, 
                                                   C=C, Bspline=Bspline, 
                                                   data_mat = data_mat, 
                                                   theta = pert_l, 
                                                   gamma = pl_params$gamma, 
                                                   p =  pl_params$p_at_conv)
          I_theta[k,l] <- I_theta[l,k] <- od_loglik_pert_both - od_loglik_pert_k - od_loglik_pert_l + od_loglik_theta # symmetry of covariance 
        }
      }
    }
    I_theta <- h_n^(-2) * I_theta
    cov_theta <- -solve(I_theta)
    # ------------------------- Estimate Cov(theta) using profile likelihood
    
    return(list(Coefficients = data.frame(Coefficient = new_theta, 
                                          SE = sqrt(diag(cov_theta))),
                converged = CONVERGED,
                converged_msg = ""))
  }
}