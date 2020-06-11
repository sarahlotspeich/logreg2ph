pl_theta_Yonly <- function(k, theta, h_n, n, N, Y_unval, Y_val, X_val, X_val_interact, C, comp_dat_all, gamma0)
{
  pert <- theta; pert[k] <- pert[k] + h_n
  pl_params <- profile_out_Yonly(theta = pert, n = n, N = N, Y_unval=Y_unval, Y_val=Y_val, X_val=X_val, X_val_interact = X_val_interact,
                                 C=C, comp_dat_all = comp_dat_all, gamma0 = gamma0)
  if(pl_params$converged)
  {
    od_loglik_pert <- observed_data_loglik(N = N, n = n, Y_unval=Y_unval, Y_val=Y_val, X_val=X_val, X_val_interact = X_val_interact, C=C, 
                                           comp_dat_all = comp_dat_all, theta = pert, gamma = pl_params$gamma)
    
    return(od_loglik_pert)
  }
}