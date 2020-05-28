VERBOSE <- FALSE
n <- 1000 #500; 750; 1000; 2500; 5000
pv <- 0.25
nv <- n * pv
prevZ <- 0.15
beta0 <- -0.5; beta1 <- -1.5182; beta2 <- 0.9428
theta0 <- -1.3241; theta1 <- 0.4361; theta2 <- 3.3096; theta3 <- -1.0773; theta4 <- 0.1268 

Y_unval = "Ystar"; X_unval = "binXstar"; C = NULL; Validated = "V"
Y_val = "Y"; X_val = "binX"; Bspline = paste0("bs", seq(1, 2)); MAX_ITER = 5000; TOL_theta = 1E-6; TOL_gamma = 1E-4; TOL_p = 1E-4

sim_seed <- sample(x = seq(1,1000000), size = 1)
sim_res <- data.frame(rep = integer(), beta0 = numeric(), beta0_se = numeric(), beta1 = numeric(), beta1_se = numeric())

for (s in 1:100)
{
  X <- rnorm(n = n, mean = 15, sd = 6)
  binX <- as.numeric(X >= sqrt(200))
  Y <- rbinom(n = n, size = 1, prob = (1 + exp(-(beta0 + beta1*binX)))^(-1))
  Xstar <- rnorm(n = n, mean = X, sd = sqrt(15))
  binXstar <- as.numeric(Xstar >= sqrt(200))
  Ystar <- rbinom(n = n, size = 1, prob = (1 + exp(-(theta0 + theta1*binXstar + theta2*Y + theta3*binX)))^(-1))
  
  B <- matrix(NA, nrow = n, ncol = 2)
  B[,1] <- as.numeric(binXstar == 0)
  B[,2] <- as.numeric(binXstar == 1)
  colnames(B) <- paste0("bs", seq(1, 2))
  
  # Y <- rbinom(n = n, size = 1, 
  #             prob = (1 + exp(-(beta0 + beta1*binX + beta2*Z)))^(-1))
  # Ystar <- rbinom(n = n, size = 1, 
  #                 prob = (1 + exp(-(theta0 + theta1*binXstar + theta2*Y + theta3*binX + theta4*Z)))^(-1))
  # B <- matrix(NA, nrow = n, ncol = 2*2)
  # B[,1] <- as.numeric(binXstar == 0 & Z == 0)
  # B[,2] <- as.numeric(binXstar == 0 & Z == 1)
  # B[,3] <- as.numeric(binXstar == 1 & Z == 0)
  # B[,4] <- as.numeric(binXstar == 1 & Z == 1)
  # colnames(B) <- paste0("bs", seq(1, 2*2))
  
  V <- seq(1,n) %in% sample(seq(1,n), nv, replace = FALSE)
  data <- data.frame(Y, binX, Ystar, binXstar, Z, B, V)
  
  n <- nrow(data)
  n_v <- sum(data[,Validated])
  
  # Reorder so that the n_v validated subjects are first ------------
  data <- data[order(as.numeric(data[,Validated]), decreasing = TRUE),]
  
  # Add the B spline basis ------------------------------------------
  sn <- ncol(data[,Bspline])
  if(0 %in% colSums(data[c(1:n_v),Bspline]))
  {
    return(list(Coefficients = data.frame(Coefficient = NA, 
                                          SE = NA),
                converged = FALSE,
                converged_msg = "Empty sieve in validated data."))
  }
  # ------------------------------------------ Add the B spline basis
  
  # Save distinct X -------------------------------------------------
  x_obs <- data.frame(unique(data[1:n_v,c(X_val)]))
  x_obs <- data.frame(x_obs[order(x_obs[,1]),])
  colnames(x_obs) <- c(X_val)
  m <- nrow(x_obs)
  x_obs_stacked <- data.frame(x_obs[rep(seq(1,m), each = (n-n_v)),])
  # x_obs_stacked <- do.call(rbind, replicate(n = (n-n_v), expr = x_obs, simplify = FALSE))
  # x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[,1]),])
  colnames(x_obs_stacked) <- c(X_val)
  #suppressMessages(data %<>% dplyr::left_join(data.frame(x_k, k = 1:m)))
  
  # Save static (X*,Y*,X,Y,C) since they don't change ---------------
  comp_dat_val <- data[c(1:n_v),c(Y_unval, X_unval, C, Bspline, X_val, Y_val)]
  comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
  comp_dat_val <- comp_dat_val[,c(Y_unval, X_unval, C, Bspline, X_val, Y_val, "k")]
  comp_dat_val <- data.matrix(comp_dat_val)
  # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
  # one row per x) --------------------------------------------------
  suppressWarnings(comp_dat_unval <- cbind(data[-c(1:n_v),c(Y_unval, X_unval, C, Bspline)],
                                           x_obs_stacked))
  comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
  comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
  colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
  comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1), 
                                      k = rep(rep(seq(1,m), each = (n-n_v)), times = 2)))
  
  comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  
  num_pred <- length(X_val) + length(C) #preds in analysis model --
  prev_theta <- theta0 <- matrix(0, nrow = (num_pred+1), ncol = 1)
  prev_gamma <- gamma0 <- matrix(0, nrow = (length(Y_val) + length(X_unval) + num_pred + 1), ncol = 1)
  
  theta_formula <- as.formula(paste0(Y_val, "~", paste(c(X_val, C), collapse = "+")))
  theta_design_mat <- cbind(int = 1, comp_dat_all[,c(X_val,C)])
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[,c(X_unval, Y_val, X_val, C)])
  # If unvalidated variable was left blank, assume error-free -------
  ## Need to write simplification here 
  
  # Initialize B-spline coefficients {p_kj}  ------------ 
  ## Numerators sum B(Xi*) over k = 1,...,m -------------
  ## Save as p_val_num for updates ----------------------
  ## (contributions don't change) -----------------------
  p_val_num <- rowsum(x = comp_dat_val[,Bspline], group = comp_dat_val[,"k"], reorder = TRUE)
  prev_p <- p0 <-  t(t(p_val_num)/colSums(p_val_num))
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1
  
  save_theta <- as.vector(theta0)
  
  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    pY_X <- 1/(1 + exp(-as.numeric(theta_design_mat[-c(1:n_v),] %*% prev_theta)))
    pY_X[which(comp_dat_unval[,Y_val] == 0)] <- 1-pY_X[which(comp_dat_unval[,Y_val] == 0)]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[-c(1:n_v),] %*% prev_gamma)))
    pYstar[which(comp_dat_unval[,Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_unval[,Y_unval] == 0)]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    pX <- do.call(rbind, replicate(n = (n-n_v), expr = prev_p, simplify = FALSE))
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    pX <- pX[order(rep(seq(1,m), times = (n-n_v))),]
    pX <- rbind(pX, pX)
    ### ---------------------------------------------------------- p_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)p_kjB(X*) ------------------------------
    psi_t_num <- pY_X * pYstar * comp_dat_unval[,Bspline] * pX
    ### Update denominator --------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ------------------
    #### reorder = TRUE returns them in ascending order of i ----------
    #### (rather than in order of encounter) --------------------------
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = rep(seq(1,(n-n_v)), times = 2*m))
    #### Then sum over the sn splines ---------------------------------
    #### Same ordering as psi_t_num_sumover_k_y, just only 1 column ---
    psi_t_num_sumover_k_y_j <- rowSums(psi_t_num_sumover_k_y)
    #### Replicate psi_t_num_sumover_k_y_j as the denominator ---------
    psi_t_denom <- matrix(rep(rep(psi_t_num_sumover_k_y_j, times = 2*m), sn), nrow = nrow(psi_t_num), ncol = ncol(psi_t_num), byrow = FALSE)
    ### And divide them! ----------------------------------------------
    psi_t <- psi_t_num/psi_t_denom
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
    ## Update the w_kyi for unvalidated subjects ----------------------
    ## by summing across the splines/ columns of psi_t ----------------
    ## w_t is ordered by i = 1, ..., n --------------------------------
    w_t <- rowSums(psi_t)
    ## ---------------------- Update the w_kyi for unvalidated subjects
    ## For validated subjects, w_t = I(Xi=xk) so make them all 0 ------
    #w_t[rep(data[,Validated], 2*m)] <- 0
    ## then place a 1 in the w_t_val positions ------------------------
    #w_t[w_t_val] <- 1
    ## Check: w_t sums to 1 over within i -----------------------------
    # table(rowSums(rowsum(w_t, group = rep(rep(seq(1,n), times = m), times = 2))))
    
    ## Update the u_kji for unvalidated subjects ----------------------
    ## by summing over Y = 0/1 w/i each i, k --------------------------
    ## add top half of psi_t (y = 0) to bottom half (y = 1) -----------
    u_t <- psi_t[c(1:(m*(n-n_v))),] + psi_t[-c(1:(m*(n-n_v))),]
    ## make u_t for the (1:n_v) validated subjects = 0 ----------------
    ## so that they won't contribute to updated p_kj ------------------
    #u_t[rep(data[,Validated], times = m),] <- 0
    ## ---------------------- Update the u_kji for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted logistic regression ----------------
    ### Gradient ------------------------------------------------------
    mu <- theta_design_mat %*% prev_theta
    w_t <- c(rep(1,n_v), w_t)
    gradient_theta <- matrix(data = c(colSums(w_t * c((comp_dat_all[,Y_val]-1 + exp(-mu)/(1+exp(-mu)))) * theta_design_mat)), ncol = 1)
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
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    if(VERBOSE) print(new_theta)
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta)<TOL_theta
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    mu <- gamma_design_mat %*% prev_gamma
    gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
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
      suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    #if(VERBOSE) print(new_gamma)
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma)<TOL_gamma
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., n -----------
    new_p_num <- p_val_num + 
      rowsum(u_t, group = rep(seq(1,m), each = (n-n_v)), reorder = TRUE)
    new_p <- t(t(new_p_num)/colSums(new_p_num))
    #if(VERBOSE) print(new_p[1,])
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p)<TOL_p
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(VERBOSE & it%%25 == 0) print(paste("Iteration", it, "complete."))
    
    save_theta %<>% append(new_theta)
    
    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p 
  }
  
  # observed-data log-likelihood --------------------------------------
  mu_theta <- theta_design_mat %*% new_theta
  pY_X <- (1 + exp(-mu_theta))^(-1)
  pY_X[which(comp_dat_all[,Y_val] == 0)] <- 1 - pY_X[which(comp_dat_all[,Y_val] == 0)] 
  unique(data.frame(comp_dat_all[,c(Y_val, X_val)], pY_X)) %>% 
    dplyr::group_by(binX) %>% dplyr::summarise(sum_pY_X = sum(pY_X))
  
  mu_gamma <- gamma_design_mat %*% new_gamma
  pYstar<- (1 + exp(-mu_gamma))^(-1)
  pYstar[which(comp_dat_all[,Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[,Y_unval] == 0)] 
  unique(data.frame(comp_dat_all[,c(Y_val, X_val, Y_unval, X_unval)], pYstar)) %>% 
    dplyr::group_by(binXstar, Y, binX) %>% dplyr::summarise(sum_pYstar = sum(pYstar))
  
  p_k <- new_p[comp_dat_all[,"k"],]
  
  # Validated subjects
  obs_loglik <- sum(log(pY_X[c(1:n_v)])) + sum(log(pYstar[c(1:n_v)])) + sum(log(p_k[c(1:n_v),])*comp_dat_all[c(1:n_v),Bspline])
  
  unval_summand <- rowSums(rowsum(x = c(pY_X[-c(1:n_v)]*pYstar[-c(1:n_v)])*p_k[-c(1:n_v),]*comp_dat_all[-c(1:n_v),Bspline], group = rep(seq(1,(n-n_v)), times = 2*m)))
  obs_loglik <- obs_loglik + sum(log(unval_summand))
  
  od_loglik_at_SMLE <- observed_data_loglik(n = n, n_v = n_v, 
                                            Y_unval=Y_unval, Y_val=Y_val, 
                                            X_unval=X_unval, X_val=X_val, 
                                            C=C, Bspline=Bspline, 
                                            comp_dat_all = comp_dat_all, 
                                            theta = new_theta,
                                            gamma = new_gamma, 
                                            p = new_p)
  
  od_ll <- vector()
  
  for (t in seq(1,length(save_theta),by=2))
  {
    od_ll %<>% append(observed_data_loglik(n = n, n_v = n_v, 
                                           Y_unval=Y_unval, Y_val=Y_val, 
                                           X_unval=X_unval, X_val=X_val, 
                                           C=C, Bspline=Bspline, 
                                           comp_dat_all = comp_dat_all, 
                                           theta = save_theta[t:(t+1)],
                                           gamma = new_gamma, 
                                           p = new_p))
  }
  
  data.frame(it = 1:length(od_ll), loglik = od_ll) %>% ggplot(aes(x = it, y = loglik)) + geom_line()
  
  pl_theta <- observed_data_loglik(n = n, n_v = n_v, 
                                   Y_unval=Y_unval, Y_val=Y_val, 
                                   X_unval=X_unval, X_val=X_val, 
                                   C=C, Bspline=Bspline, 
                                   comp_dat_all = comp_dat_all, 
                                   theta = new_theta,
                                   gamma = new_gamma, 
                                   p = new_p)
  
  h_n <- n^(-1/2) # perturbation ----------------------------
  I_theta <- matrix(0, nrow = nrow(new_theta), ncol = nrow(new_theta))
  
  pert_alpha <- new_theta; pert_alpha[1] <- pert_alpha[1] + h_n
  pert_alpha2 <- pert_alpha; pert_alpha2[1] <- pert_alpha2[1] + h_n
  pert_beta <- new_theta; pert_beta[2] <- pert_beta[2] + h_n
  pert_beta2 <- pert_beta; pert_beta2[2] <- pert_beta2[2] + h_n
  pert_both <- new_theta + h_n
  
  pl_params <- profile_out(theta = pert_alpha, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num,
                           newton_step_scale = newton_step_scale)
  pl_pert_alpha <- observed_data_loglik(n = n, n_v = n_v, 
                                        Y_unval=Y_unval, Y_val=Y_val, 
                                        X_unval=X_unval, X_val=X_val, 
                                        C=C, Bspline=Bspline, 
                                        comp_dat_all = comp_dat_all, 
                                        theta = pert_alpha,
                                        gamma = pl_params$gamma_at_conv, 
                                        p = pl_params$p_at_conv)
  
  pl_params <- profile_out(theta = pert_alpha2, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num,
                           newton_step_scale = newton_step_scale)
  pl_pert_alpha2 <- observed_data_loglik(n = n, n_v = n_v, 
                                         Y_unval=Y_unval, Y_val=Y_val, 
                                         X_unval=X_unval, X_val=X_val, 
                                         C=C, Bspline=Bspline, 
                                         comp_dat_all = comp_dat_all, 
                                         theta = pert_alpha2,
                                         gamma = pl_params$gamma_at_conv, 
                                         p = pl_params$p_at_conv)
  
  pl_params <- profile_out(theta = pert_beta, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num,
                           newton_step_scale = newton_step_scale)
  pl_pert_beta <- observed_data_loglik(n = n, n_v = n_v, 
                                       Y_unval=Y_unval, Y_val=Y_val, 
                                       X_unval=X_unval, X_val=X_val, 
                                       C=C, Bspline=Bspline, 
                                       comp_dat_all = comp_dat_all, 
                                       theta = pert_beta,
                                       gamma = pl_params$gamma_at_conv, 
                                       p = pl_params$p_at_conv)
  
  pl_params <- profile_out(theta = pert_beta2, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num,
                           newton_step_scale = newton_step_scale)
  pl_pert_beta2 <- observed_data_loglik(n = n, n_v = n_v, 
                                        Y_unval=Y_unval, Y_val=Y_val, 
                                        X_unval=X_unval, X_val=X_val, 
                                        C=C, Bspline=Bspline, 
                                        comp_dat_all = comp_dat_all, 
                                        theta = pert_beta2,
                                        gamma = pl_params$gamma_at_conv, 
                                        p = pl_params$p_at_conv)
  
  pl_params <- profile_out(theta = pert_both, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num,
                           newton_step_scale = newton_step_scale)
  pl_pert_both <- observed_data_loglik(n = n, n_v = n_v, 
                                       Y_unval=Y_unval, Y_val=Y_val, 
                                       X_unval=X_unval, X_val=X_val, 
                                       C=C, Bspline=Bspline, 
                                       comp_dat_all = comp_dat_all, 
                                       theta = pert_both,
                                       gamma = pl_params$gamma_at_conv, 
                                       p = pl_params$p_at_conv)
  
  I_theta[1,1] <- pl_pert_alpha2 - 2*pl_pert_alpha + pl_theta
  I_theta[2,2] <- pl_pert_beta2 - 2*pl_pert_beta + pl_theta
  I_theta[1,2] <- I_theta[2,1] <- pl_pert_both - pl_pert_alpha - pl_pert_beta + pl_theta
  
  I_theta <- (1/(h_n^2))*I_theta
  cov_theta <- -solve(I_theta)
  
  sim_res <- rbind(sim_res, data.frame(rep = s, beta0 = new_theta[1], beta0_se = sqrt(diag(cov_theta))[1], beta1 = new_theta[2], beta1_se = sqrt(diag(cov_theta))[2]))
  write.csv(sim_res, paste0("seed", sim_seed, ".csv"), row.names = FALSE)
}

