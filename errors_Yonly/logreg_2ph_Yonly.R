logreg_2ph_Yonly <- function(N, n, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, X_val_interact = NULL,
                              data, initial_lr_params, h_n_scale, noSE, TOL, MAX_ITER, rescale)
{
  # Save static (Y*,X,Y,C) since they don't change ---------------
  comp_dat_val <- data[c(1:n),c(Y_unval, C, X_val, X_val_interact, Y_val)]
  comp_dat_val <- comp_dat_val[,c(Y_unval, C, X_val, X_val_interact, Y_val)]
  comp_dat_val <- data.matrix(comp_dat_val)
  
  # 2 Nxd matrices (y=0/y=1) of each (one column per person, --
  # one row per x) --------------------------------------------------
  comp_dat_unval <- re_data[-c(1:n),c(Y_unval, C, X_val, X_val_interact)]
  comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
  comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
  colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
  comp_dat_unval <- data.matrix(rbind(comp_dat_y0, comp_dat_y1))
  
  comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  
  theta_formula <- as.formula(paste0(Y_val, "~", paste(c(X_val, gsub("x", "*", X_val_interact), C), collapse = "+")))
  theta_design_mat <- cbind(int = 1, comp_dat_all[,c(X_val,X_val_interact,C)])
  
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[,c(Y_val, X_val, C)])
  
  # Standardize Y_val in gamma_design_mat to N(0,1) -----------------
  ## Everything else was already standardized above -----------------
  if(rescale)
  {
    gamma_design_mat[,Y_val] <- (gamma_design_mat[,Y_val] - mean(gamma_design_mat[,Y_val]))/sd(gamma_design_mat[,Y_val])  
  }
  # ----------------- Standardize Y_val in gamma_design_mat to N(0,1)
  
  # Initialize parameter values theta, gamma ------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data")))
  {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params <- "Zero"
  }
  if(initial_lr_params == "Zero")
  {
    num_pred <- length(X_val) + length(X_val_interact) + length(C) 
    prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
    prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)
  }
  if(initial_lr_params == "Complete-data")
  {
    prev_theta <- theta0 <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(re_data[c(1:n),]))$coefficients, ncol = 1)
    prev_gamma <- gamma0 <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(re_data[c(1:n),]))$coefficient, ncol = 1)
  }

  y0 <- comp_dat_unval[,Y_val] == 0
  ystar0 <- comp_dat_unval[,Y_unval] == 0
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1
  
  # Estimate using EM -------------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi for unvalidated subjects ------------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric(theta_design_mat[-c(1:n),] %*% prev_theta)
    pY_X <- 1/(1 + exp(-mu_theta)) 
    pY_X[y0] <- 1-pY_X[y0]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[-c(1:n),] %*% prev_gamma)))
    pYstar[ystar0] <- 1 - pYstar[ystar0]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|Y,X,C) ------------------------------------------
    psi_t_num <- pY_X * pYstar
    ### Update denominator --------------------------------------------
    psi_t_denom <- rowsum(psi_t_num, group = rep(seq(1,(N-n)), times = 2))
    #### Avoid NaN resulting from dividing by 0 -----------------------
    psi_t_denom[psi_t_denom == 0] <- 1
    ### And divide them! ----------------------------------------------
    w_t <- psi_t_num/psi_t_denom
    ## ------------------------ Update the psi for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------
    
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted logistic regression ----------------
    ### Gradient ------------------------------------------------------
    mu <- theta_design_mat %*% prev_theta
    w_t <- c(rep(1,n), w_t)
    gradient_theta <- matrix(data = c(colSums(w_t * c((comp_dat_all[,Y_val]-1 + exp(-mu)/(1+exp(-mu)))) * theta_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c((exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))*w_t*theta_design_mat
    hessian_theta <- apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    
    new_theta <- tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta,
                          error = function(err) {matrix(NA, nrow = nrow(prev_theta))
                          })
    if (TRUE %in% is.na(new_theta))
    {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta)<TOL
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    mu <- gamma_design_mat %*% prev_gamma
    gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[,c(Y_unval)]-1 + exp(-mu)/(1+exp(-mu)))) * gamma_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c(w_t*(exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))*gamma_design_mat
    hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
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
    # M Step ----------------------------------------------------------
    
    all_conv <- c(theta_conv, gamma_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_gamma <- new_gamma
  }
  # ------------------------------------------------- Estimate using EM
  
  if(!CONVERGED) 
  {
    if(it > MAX_ITER)
    {
      CONVERGED_MSG = "MAX_ITER reached"
    }
    return(list(Coefficients = data.frame(Coefficient = matrix(NA, nrow = nrow(prev_theta)), 
                                          SE = NA),
                h_n = NA,
                converged = CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = NA,
                initial_vals = initial_lr_params))
  }
  if(CONVERGED)
  {
    CONVERGED_MSG <- "Converged"
  }
  
  if(noSE)
  {
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik_Yonly(N = N, n = n, Y_unval=Y_unval, Y_val=Y_val, X_val=X_val,
                                                  X_val_interact = X_val_interact, C=C, comp_dat_all = comp_dat_all,
                                                  theta = new_theta, gamma = new_gamma)
    
    rownames(new_theta) <- c("Intercept", colnames(theta_design_mat)[-1])
    re_theta <- new_theta
    rownames(new_gamma) <- c("Intercept", X_unval, Y_val, X_val, C)
    re_gamma <- new_gamma
    if(rescale)
    {
      re_theta[c(2:(1+length(c(X_val,X_val_interact,C))))] <- re_theta[c(2:(1+length(c(X_val,X_val_interact,C))))]/re_scale[c(2:(1+length(c(X_val,X_val_interact,C))))]
      re_theta[1] <- re_theta[1] - sum(re_theta[c(2:(1+length(c(X_val,X_val_interact,C))))]*re_shift[c(2:(1+length(c(X_val,X_val_interact,C))))])
      
      re_gamma[c(2:(1+length(c(X_val,C))))] <- re_gamma[c(2:(1+length(c(X_val,C))))]/re_scale[c(2:(1+length(c(X_val,C))))]
      re_gamma[1] <- re_gamma[1] - sum(re_gamma[c(2:(1+length(c(X_val,C))))]*re_shift[c(2:(1+length(c(X_val,C))))])
    }
    
    return(list(Coefficients = data.frame(Coefficient = re_theta, 
                                          SE = NA),
                h_n = NA,
                converged = CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                gamma_at_conv = re_gamma, 
                od_loglik_at_conv = od_loglik_theta,
                initial_vals = initial_lr_params))
  } else
  {
    rownames(new_theta) <- c("Intercept", colnames(theta_design_mat)[-1])
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_n <- h_n_scale*N^(-1/2) # perturbation ----------------------------
    
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik_Yonly(N = N, n = n, Y_unval = Y_unval, Y_val = Y_val, X_val = X_val, X_val_interact = X_val_interact, C = C, 
                                            comp_dat_all = comp_dat_all, theta = new_theta, gamma = new_gamma)
    
    I_theta <- matrix(od_loglik_theta, nrow = nrow(new_theta), ncol = nrow(new_theta))
    
    #single perts: 271.658 sec elapsed
    single_pert_theta <- sapply(X = seq(1,ncol(I_theta)), FUN = pl_theta_Yonly, theta = new_theta, h_n = h_n, n = n, N = N, 
                                Y_unval = Y_unval, Y_val = Y_val, X_val = X_val, X_val_interact = X_val_interact, C = C, 
                                comp_dat_all = comp_dat_all, gamma0 = new_gamma)
    
    spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)), ncol = ncol(I_theta), byrow = FALSE)
    #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
    I_theta <- I_theta - spt_wide - t(spt_wide)
    
    for (c in 1:ncol(I_theta))
    {
      pert_theta <- new_theta; pert_theta[c] <- pert_theta[c] + h_n
      ## Avoid 
      double_pert_theta <- sapply(X = seq(c,ncol(I_theta)), FUN = pl_theta_Yonly, theta = pert_theta, h_n = h_n, n = n, N = N, 
                                  Y_unval=Y_unval, Y_val=Y_val, X_val=X_val, X_val_interact = X_val_interact, C=C,  
                                  comp_dat_all = comp_dat_all, gamma0 = new_gamma)
      dpt <- matrix(0, nrow = nrow(I_theta), ncol = ncol(I_theta))
      dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal 
      if(c < ncol(I_theta))
      {
        dpt[c,-(1:c)] <- dpt[-(1:c),c] <- double_pert_theta[-1] #And fill the others in on the cth row/ column 
      }
      I_theta <- I_theta + dpt
    }
    
    I_theta <- h_n^(-2) * I_theta
    
    tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
             error = function(err) {matrix(NA, nrow = nrow(prev_gamma))
             })
    cov_theta <- -solve(I_theta)
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # Scale everything back ------------------------------------------------
    re_theta <- new_theta
    if(TRUE %in% (diag(cov_theta) < 0))
    {
      warning("Negative variance estimate. Increase the h_n_scale parameter and repeat variance estimation.")
      SE_cov <- FALSE
    }
    re_se_theta <- tryCatch(expr = sqrt(diag(cov_theta)), 
                            warning = function(w) {matrix(NA, nrow = nrow(prev_theta))}) 
    re_gamma <- new_gamma
    
    if(rescale)
    {
      re_theta[c(2:(1+length(c(X_val,X_val_interact,C))))] <- re_theta[c(2:(1+length(c(X_val,X_val_interact,C))))]/re_scale[c(2:(1+length(c(X_val,X_val_interact,C))))]
      re_theta[1] <- re_theta[1] - sum(re_theta[c(2:(1+length(c(X_val,X_val_interact,C))))]*re_shift[c(2:(1+length(c(X_val,X_val_interact,C))))])
      
      re_se_theta[c(2:(1+length(c(X_val,X_val_interact,C))))] <- re_se_theta[c(2:(1+length(c(X_val,X_val_interact,C))))]/re_scale[c(2:(1+length(c(X_val,X_val_interact,C))))]
      re_se_theta[1] <- cov_theta[1,1] + sum(diag(cov_theta)[c(2:(1+length(c(X_val,X_val_interact,C))))]*(re_shift[c(2:(1+length(c(X_val,X_val_interact,C))))]/re_scale[c(2:(1+length(c(X_val,X_val_interact,C))))])^2)
      for (p1 in 1:ncol(cov_theta))
      {
        for (p2 in p1:ncol(cov_theta))
        {
          if(p1 < p2 & p1 == 1)
          {
            re_se_theta[1] <- re_se_theta[1] - (re_shift[p2]/re_scale[p2])*cov_theta[p1,p2]
          }
          if(p1 < p2 & p1 > 1)
          {
            re_se_theta[1] <- re_se_theta[1] + (re_shift[p1]/re_scale[p1])*(re_shift[p2]/re_scale[p2])*cov_theta[p1,p2]
          }
        }
      }
      re_se_theta[1] <- sqrt(re_se_theta[1])
      
      re_gamma[c(2:(1+length(c(X_val,C))))] <- re_gamma[c(2:(1+length(c(X_val,C))))]/re_scale[c(2:(1+length(c(X_val,C))))]
      re_gamma[1] <- re_gamma[1] - sum(re_gamma[c(2:(1+length(c(X_val,C))))]*re_shift[c(2:(1+length(c(X_val,C))))])
    }
    # ------------------------------------------------ Scale everything back
    return(list(Coefficients = data.frame(Coefficient = re_theta, 
                                          SE = re_se_theta),
                h_n = h_n,
                I_theta = I_theta,
                converged = CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                gamma_at_conv = re_gamma, 
                od_loglik_at_conv = od_loglik_theta,
                initial_vals = initial_lr_params))
  }
}