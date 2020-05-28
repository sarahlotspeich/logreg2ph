profile_out_all <- function(theta, n_v, n, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, Bspline=NULL, comp_dat_all, gamma0, p0, p_val_num, newton_step_scale = 1, TOL = 1E-4, MAX_ITER = 5000)
{
  sn <- ncol(p0)
  m <- nrow(p0)
  
  prev_theta <- theta
  prev_gamma <- gamma0
  prev_p <- p0
  
  theta_design_mat <- cbind(int = 1, comp_dat_all[,c(X_val,C)])
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[,c(X_unval, Y_val, X_val, C)])
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1
  # Estimate gamma/p using EM ----------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    pY_X <- 1/(1 + exp(-as.numeric(theta_design_mat[-c(1:n_v),] %*% prev_theta)))
    pY_X[which(comp_dat_all[-c(1:n_v),Y_val] == 0)] <- 1-pY_X[which(comp_dat_all[-c(1:n_v),Y_val] == 0)]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[-c(1:n_v),] %*% prev_gamma)))
    pYstar[which(comp_dat_all[-c(1:n_v),Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[-c(1:n_v),Y_unval] == 0)]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    pX <- do.call(rbind, replicate(n = (n-n_v), expr = prev_p, simplify = FALSE))
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    pX <- pX[order(rep(seq(1,m), times = (n-n_v))),]
    pX <- rbind(pX, pX)
    ### ---------------------------------------------------------- p_kj
    ###################################################################
    ### Update the psi_kyji for unvalidated subjects ------------------
    psi_t_num <- pY_X * pYstar * comp_dat_all[-c(1:n_v),Bspline] * pX
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = rep(seq(1,(n-n_v)), times = 2*m))
    psi_t_num_sumover_k_y_j <- rowSums(psi_t_num_sumover_k_y)
    psi_t_denom <- matrix(rep(rep(psi_t_num_sumover_k_y_j, times = 2*m), sn), nrow = nrow(psi_t_num), ncol = ncol(psi_t_num), byrow = FALSE)
    psi_t <- psi_t_num/psi_t_denom
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
    ## Update the w_kyi for unvalidated subjects ----------------------
    w_t <- rowSums(psi_t)
    ## ---------------------- Update the w_kyi for unvalidated subjects
    ## Update the u_kji for unvalidated subjects ---------------------- 
    u_t <- psi_t[c(1:(m*(n-n_v))),] + psi_t[-c(1:(m*(n-n_v))),]
    ## ---------------------- Update the u_kji for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
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
    new_theta[2] <- theta[2]
    if (TRUE %in% is.na(new_theta))
    {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      new_theta[2] <- theta[2]
    }
    theta_conv <- abs(new_theta - prev_theta)<TOL
    
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
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma)<TOL
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., n -----------
    new_p_num <- p_val_num + 
      rowsum(u_t, group = rep(seq(1,m), each = (n-n_v)), reorder = TRUE)
    new_p <- t(t(new_p_num)/colSums(new_p_num))
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p)<TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p 
  }
  if(it == MAX_ITER & !CONVERGED) CONVERGED_MSG <- "MAX_ITER reached"
  if(CONVERGED) CONVERGED_MSG <- "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("theta_at_conv" = new_theta,
              "gamma_at_conv" = new_gamma,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
