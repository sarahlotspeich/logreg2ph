observed_data_loglik <- function(n, n_v, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, comp_dat_all, theta, gamma, p)
{
  sn <- ncol(p)
  m <- nrow(p)
  
  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n_v),c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(comp_dat_all[c(1:n_v),c(Y_val)]) == 0, 1-pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over log[P(Yi*|Xi*,Yi,Xi)] -----------------------------------------------
  pYstar <- 1/(1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n_v),c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar <- ifelse(as.vector(comp_dat_all[c(1:n_v),Y_unval]) == 0, 1-pYstar, pYstar)
  return_loglik <- return_loglik + sum(log(pYstar))
  ## ----------------------------------------------- Sum over log[P(Yi*|Xi*,Yi,Xi)]
  #################################################################################
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
  pX <- p[comp_dat_all[c(1:n_v),"k"],]
  log_pX <- log(pX)
  log_pX[log_pX == -Inf] <- 0
  return_loglik <- return_loglik + sum(comp_dat_all[c(1:n_v),Bspline] * log_pX)
  ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  #################################################################################
  # -------------------------------------------------------- For validated subjects
  
  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n_v),c(X_val, C)]) %*% theta)))
  pY_X[which(comp_dat_all[-c(1:n_v),Y_val] == 0)] <- 1-pY_X[which(comp_dat_all[-c(1:n_v),Y_val] == 0)]
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate P(Yi*|Xi*,y,xk) for all (y,xk) ------------------------------------
  pYstar <- 1/(1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n_v),c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar[which(comp_dat_all[-c(1:n_v),Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[-c(1:n_v),Y_unval] == 0)]
  ## ------------------------------------ Calculate P(Yi*|Xi*,y,xk) for all (y,xk)
  ################################################################################
  ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
  pX <- p[comp_dat_all[-c(1:n_v),"k"],]
  ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x p_kj ---------------------------
  person_sum <- rowsum(c(pY_X*pYstar)*comp_dat_all[-c(1:n_v),Bspline]*pX, group = rep(seq(1,(n-n_v)), times = 2*m))
  person_sum <- rowSums(person_sum)
  log_person_sum <- log(person_sum)
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}