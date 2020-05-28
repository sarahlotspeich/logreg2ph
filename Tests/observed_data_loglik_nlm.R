observed_data_loglik_nlm <- function(params, n, n_v, m, sn, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, comp_dat_all)
{
  p_theta <- length(c(X_val, C))+1
  theta <- matrix(params[c(1:p_theta)], ncol = 1)
  p_gamma <- length(c(X_unval, Y_val, X_val, C))+1
  gamma <- matrix(params[p_theta+c(1:p_gamma)], ncol = 1)
  p_num <- matrix(params[-c(1:(p_theta+p_gamma))], nrow = m, ncol = sn)
  p <- t(t(p_num)/colSums(p_num))

  mu_theta <- cbind(1, comp_dat_all[,c(X_val, C)]) %*% theta
  pY_X <- (1 + exp(-mu_theta))^(-1)
  pY_X[which(comp_dat_all[,Y_val] == 0)] <- 1 - pY_X[which(comp_dat_all[,Y_val] == 0)] 
  
  mu_gamma <- cbind(1, comp_dat_all[,c(X_unval, Y_val, X_val, C)]) %*% gamma
  pYstar<- (1 + exp(-mu_gamma))^(-1)
  pYstar[which(comp_dat_all[,Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[,Y_unval] == 0)] 
  
  p_k <- p[comp_dat_all[,"k"],]
  
  # Validated subjects ----------------------------------
  obs_loglik <- sum(log(pY_X[c(1:n_v)])) + sum(log(pYstar[c(1:n_v)])) + sum(log(p_k[c(1:n_v),])*comp_dat_all[c(1:n_v),Bspline])
  
  # Unvalidated subjects --------------------------------
  unval_summand <- rowSums(rowsum(x = c(pY_X[-c(1:n_v)]*pYstar[-c(1:n_v)])*p_k[-c(1:n_v),]*comp_dat_all[-c(1:n_v),Bspline], group = rep(seq(1,(n-n_v)), times = 2*m)))
  obs_loglik <- obs_loglik + sum(log(unval_summand))
  return(-obs_loglik)
}
