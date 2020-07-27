pl_theta_Yonly <- function(k, theta, h_n, n, N, Y_unval, Y_val, X_val, X_val_interact, C, comp_dat_all, gamma0, TOL, MAX_ITER)
{
  pert <- theta; pert[k] <- pert[k] + h_n
  pl_params <- profile_out_Yonly(theta = pert, n = n, N = N, Y_unval = Y_unval, Y_val = Y_val, 
                                 X_val = X_val, X_val_interact = X_val_interact,
                                 C = C, comp_dat_all = comp_dat_all, gamma0 = gamma0, TOL = TOL, MAX_ITER = MAX_ITER)
  if (pl_params$converged)
  {
    od_loglik_pert <- observed_data_loglik_Yonly(N = N, n = n, Y_unval = Y_unval, Y_val = Y_val, 
                                                 X_val = X_val, X_val_interact = X_val_interact, C = C, 
                                                 comp_dat_all = comp_dat_all, theta = pert, gamma = pl_params$gamma)
  } else
  {
    od_loglik_pert <- NA
  }
  return(od_loglik_pert)
}