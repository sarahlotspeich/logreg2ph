observed_data_loglik_Yonly <- function(N, n, Y_unval = NULL, Y_val = NULL, X_val = NULL, C = NULL, X_val_interact = NULL, comp_dat_all, theta, gamma) {
  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1 / (1 + exp( - as.numeric(cbind(int = 1, comp_dat_all[c(1:n), c(X_val, X_val_interact, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(comp_dat_all[c(1:n), c(Y_val)]) == 0, 1 - pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over log[P(Yi*|Xi*,Yi,Xi)] -----------------------------------------------
  pYstar <- 1 / (1 + exp( - as.numeric(cbind(int = 1, comp_dat_all[c(1:n), c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar <- ifelse(as.vector(comp_dat_all[c(1:n), Y_unval]) == 0, 1 - pYstar, pYstar)
  return_loglik <- return_loglik + sum(log(pYstar))
  ## ----------------------------------------------- Sum over log[P(Yi*|Xi*,Yi,Xi)]
  #################################################################################
  # -------------------------------------------------------- For validated subjects
  
  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1 / (1 + exp( - as.numeric(cbind(int = 1, comp_dat_all[-c(1:n), c(X_val, X_val_interact, C)]) %*% theta)))
  pY_X[which(comp_dat_all[-c(1:n), Y_val] == 0)] <- 1 - pY_X[which(comp_dat_all[-c(1:n), Y_val] == 0)]
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate P(Yi*|Xi*,y,xk) for all (y,xk) ------------------------------------
  pYstar <- 1 / (1 + exp( - as.numeric(cbind(int = 1, comp_dat_all[-c(1:n), c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar[which(comp_dat_all[-c(1:n), Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[-c(1:n), Y_unval] == 0)]
  ## ------------------------------------ Calculate P(Yi*|Xi*,y,xk) for all (y,xk)
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) -------------------------------------------
  person_sum <- rowsum(matrix(c(pY_X * pYstar), ncol = 1), group = rep(seq(1, (N-n)), times = 2))
  log_person_sum <- log(person_sum)
  log_person_sum[log_person_sum == -Inf] <- 0
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}

profile_out_Yonly <- function(theta, n, N, Y_unval = NULL, Y_val = NULL, X_val = NULL, C = NULL, X_val_interact = NULL, 
                              comp_dat_all, gamma0, TOL, MAX_ITER) {
  comp_dat_unval <- comp_dat_all[-c(1:n), ]
  
  theta_design_mat <- cbind(int = 1, comp_dat_all[, c(X_val, X_val_interact, C)])
  
  ### P(Y|X) --------------------------------------------------------
  y0 <-  comp_dat_unval[, Y_val] == 0
  mu_theta <- as.numeric(theta_design_mat[-c(1:n), ] %*% theta)
  pY_X <- 1/(1 + exp( - mu_theta)) 
  pY_X[y0] <- 1 - pY_X[y0]
  ### -------------------------------------------------------- P(Y|X)
  
  prev_gamma <- gamma0
  ystar0 <- comp_dat_unval[, Y_unval] == 0
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[, c(X_unval, Y_val, X_val, C)])
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "SE not converged"
  it <- 1
  
  # Estimate using EM -------------------------------------------------
  while (it < MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi for unvalidated subjects ------------------------
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1 / (1 + exp(-as.numeric(gamma_design_mat[-c(1:n), ] %*% prev_gamma)))
    pYstar[ystar0] <- 1 - pYstar[ystar0]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|Y,X,C) ------------------------------------------
    psi_t_num <- pY_X * pYstar
    ### Update denominator --------------------------------------------
    psi_t_denom <- rowsum(psi_t_num, group = rep(seq(1, (N-n)), times = 2))
    #### Avoid NaN resulting from dividing by 0 -----------------------
    psi_t_denom[psi_t_denom == 0] <- 1
    ### And divide them! ----------------------------------------------
    w_t <- psi_t_num / rep(psi_t_denom, times = 2)
    ## ------------------------ Update the psi for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    mu <- gamma_design_mat %*% prev_gamma
    w_t <- c(rep(1, n), w_t)
    gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp( - mu) / (1 + exp( - mu)))) * gamma_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c(w_t * (exp( - mu)/(1 + exp( - mu))) * (exp( - mu) / (1 + exp( - mu)) - 1)) * gamma_design_mat
    hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
                          error = function(err) {
                            matrix(NA, nrow = nrow(prev_gamma))
                          }
    )
    if (TRUE %in% is.na(new_gamma)) {
      suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma) < TOL
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if (mean(gamma_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_gamma <- new_gamma
  }
  # ------------------------------------------------- Estimate using EM
  
  if (it == MAX_ITER & !CONVERGED) {
    CONVERGED_MSG <- "MAX_ITER reached"
    new_gamma <- matrix(NA, nrow = nrow(new_gamma), ncol = 1)
  }
  if (CONVERGED) CONVERGED_MSG <- "converged"
  
  return(list("gamma_at_conv" = new_gamma,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}

pl_theta_Yonly <- function(k, theta, h_n, n, N, Y_unval, Y_val, X_val, X_val_interact, C, comp_dat_all, gamma0, TOL, MAX_ITER) {
  pert <- theta; pert[k] <- pert[k] + h_n
  pl_params <- profile_out_Yonly(theta = pert, n = n, N = N, Y_unval = Y_unval, Y_val = Y_val, 
                                 X_val = X_val, X_val_interact = X_val_interact,
                                 C = C, comp_dat_all = comp_dat_all, gamma0 = gamma0, TOL = TOL, MAX_ITER = MAX_ITER)
  if (pl_params$converged) {
    od_loglik_pert <- observed_data_loglik_Yonly(N = N, n = n, Y_unval = Y_unval, Y_val = Y_val, 
                                                 X_val = X_val, X_val_interact = X_val_interact, C = C, 
                                                 comp_dat_all = comp_dat_all, theta = pert, gamma = pl_params$gamma)
  } else {
    od_loglik_pert <- NA
  }
  return(od_loglik_pert)
}