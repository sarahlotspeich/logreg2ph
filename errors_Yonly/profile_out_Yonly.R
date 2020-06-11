profile_out_Yonly <- function(theta, n, N, Y_unval=NULL, Y_val=NULL, X_val=NULL, C=NULL, X_val_interact = NULL, 
                              comp_dat_all, gamma0, TOL, MAX_ITER)
{
  comp_dat_unval <- comp_dat_all[-c(1:n),]
  
  theta_design_mat <- cbind(int = 1, comp_dat_all[,c(X_val,X_val_interact,C)])
  
  ### P(Y|X) --------------------------------------------------------
  y0 <-  comp_dat_unval[,Y_val] == 0
  mu_theta <- as.numeric(theta_design_mat[-c(1:n),] %*% theta)
  pY_X <- 1/(1 + exp(-mu_theta)) 
  pY_X[y0] <- 1-pY_X[y0]
  ### -------------------------------------------------------- P(Y|X)
  
  prev_gamma <- gamma0
  ystar0 <- comp_dat_unval[,Y_unval] == 0
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[,c(X_unval, Y_val, X_val, C)])
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "SE not converged"
  it <- 1
  
  # Estimate using EM -------------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi for unvalidated subjects ------------------------
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
    
    if (mean(gamma_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_gamma <- new_gamma
  }
  # ------------------------------------------------- Estimate using EM
  
  # Estimate gamma/p using EM ----------------------------------------------
  if(it == MAX_ITER & !CONVERGED)
  {
    CONVERGED_MSG <- "MAX_ITER reached"
    new_gamma <- matrix(NA, nrow = nrow(new_gamma), ncol = 1)
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"

  return(list("gamma_at_conv" = new_gamma,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}