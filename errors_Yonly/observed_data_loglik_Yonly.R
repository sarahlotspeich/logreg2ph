observed_data_loglik_Yonly <- function(N, n, Y_unval = NULL, Y_val = NULL, X_val = NULL, C = NULL, X_val_interact = NULL, comp_dat_all, theta, gamma)
{
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