profile_out <- function(theta, n_v, n, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, X_val_interact = NULL, C=NULL, Bspline=NULL, comp_dat_all, gamma0, p0, p_val_num, TOL = 1E-4, MAX_ITER = 5000)
{
  sn <- ncol(p0)
  m <- nrow(p0)

  prev_gamma <- gamma0
  prev_p <- p0
  
  comp_dat_unval <- comp_dat_all[-c(1:n_v),]
  
  theta_design_mat <- cbind(int = 1, comp_dat_all[,c(X_val,X_val_interact,C)])
  
  ### P(Y|X) --------------------------------------------------------
  y0 <-  comp_dat_unval[,Y_val] == 0
  mu_theta <- as.numeric(theta_design_mat[-c(1:n_v),] %*% theta)
  pY_X <- 1/(1 + exp(-mu_theta)) 
  pY_X[y0] <- 1-pY_X[y0]
  ### -------------------------------------------------------- P(Y|X)

  ystar0 <- comp_dat_unval[,Y_unval] == 0
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[,c(X_unval, Y_val, X_val, C)])
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "SE not converged"
  it <- 1
  # Estimate gamma/p using EM ----------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[-c(1:n_v),] %*% prev_gamma)))
    pYstar[ystar0] <- 1 - pYstar[ystar0]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    pX <- prev_p[rep(rep(seq(1,m), each = (n-n_v)), times = 2),]
    ### ---------------------------------------------------------- p_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)p_kjB(X*) ------------------------------
    psi_t_num <- pY_X * pYstar * comp_dat_unval[,Bspline] * pX
    ### Update denominator --------------------------------------------
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = rep(seq(1,(n-n_v)), times = 2*m))
    psi_t_num_sumover_k_y_j <- rowSums(psi_t_num_sumover_k_y)
    psi_t_num_sumover_k_y_j[psi_t_num_sumover_k_y_j == 0] <- 1
    ### And divide them! ----------------------------------------------
    psi_t <- psi_t_num/psi_t_num_sumover_k_y_j
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
    ## Update gamma using weighted logistic regression ----------------
    w_t <- c(rep(1,n_v), w_t)
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
    #if(VERBOSE) print(new_gamma)
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma)<TOL
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., n -----------
    new_p_num <- p_val_num + 
      rowsum(u_t, group = rep(seq(1,m), each = (n-n_v)), reorder = TRUE)
    new_p <- t(t(new_p_num)/colSums(new_p_num))
    #if(VERBOSE) print(new_p[1,])
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p)<TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    all_conv <- c(gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    it <- it + 1
    
    # Update values for next iteration  -------------------------------
    prev_gamma <- new_gamma
    prev_p <- new_p 
  }
  if(it == MAX_ITER & !CONVERGED)
  {
    CONVERGED_MSG <- "MAX_ITER reached"
    new_gamma <- matrix(NA, nrow = nrow(new_gamma), ncol = 1)
    new_p <- matrix(NA, nrow = nrow(new_p), ncol = ncol(new_p))
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "gamma_at_conv" = new_gamma,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}

observed_data_loglik <- function(n, n_v, Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, X_val_interact = NULL, C=NULL, Bspline=NULL, comp_dat_all, theta, gamma, p)
{
  sn <- ncol(p)
  m <- nrow(p)
  
  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n_v),c(X_val, X_val_interact, C)]) %*% theta)))
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
  pY_X <- 1/(1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n_v),c(X_val, X_val_interact, C)]) %*% theta)))
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
  log_person_sum[log_person_sum == -Inf] <- 0
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}

hessian_row <- function(x, pm)
{
  return(colSums(x*pm))
}

pl_theta <- function(k, theta, h_n, n_v, n, Y_unval, Y_val, X_unval, X_val, X_val_interact, C, Bspline, 
                     comp_dat_all, gamma0, p0, p_val_num)
{
  pert <- theta; pert[k] <- pert[k] + h_n
  pl_params <- profile_out(theta = pert, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           X_val_interact = X_val_interact,
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = gamma0, p0 = p0, p_val_num = p_val_num)
  if(pl_params$converged)
  {
    od_loglik_pert <- observed_data_loglik(n = n, n_v = n_v, 
                                           Y_unval=Y_unval, Y_val=Y_val, 
                                           X_unval=X_unval, X_val=X_val, 
                                           X_val_interact = X_val_interact,
                                           C=C, Bspline=Bspline, 
                                           comp_dat_all = comp_dat_all, 
                                           theta = pert,
                                           gamma = pl_params$gamma, 
                                           p = pl_params$p_at_conv)
    
    return(od_loglik_pert)
  }
}

TwoPhase_LogReg <- function(Y_unval=NULL, Y_val=NULL, X_unval=NULL, X_val=NULL, C=NULL, interact = NULL,
                            Validated = NULL, Bspline=NULL, data, initial_lr_params = "Zero",
                            h_n_scale = NULL, noSE=FALSE, VERBOSE = FALSE, TOL = 1E-4, MAX_ITER = 1000,
                            rescale = FALSE)
{
  n <- nrow(data)
  n_v <- sum(data[,Validated])
  
  if(is.null(h_n_scale))
  {
      if(rescale)
      {
          h_n_scale <- 4
      } else
      {
          h_n_scale <- 1
      }
  }
  
  # Reorder so that the n_v validated subjects are first ------------
  data <- data[order(as.numeric(data[,Validated]), decreasing = TRUE),]
  
  # Add the B spline basis ------------------------------------------
  sn <- ncol(data[,Bspline])
  if(0 %in% colSums(data[c(1:n_v),Bspline]))
  {
    warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)
    
    return(list(Coefficients = data.frame(Coefficient = NA, 
                                          SE = NA),
                h_n = NA,
                converged = FALSE,
                converged_msg = "B-spline error",
                initial_vals = NA, 
                iterations = 0))
  }
  # ------------------------------------------ Add the B spline basis
  
  # Create interaction terms ----------------------------------------
  X_val_interact <- NULL
  if(length(interact) > 0)
  {
    for (int in 1:length(interact))
    {
      interact_term <- data[,interact[[int]][1]]
      for (t in 2:length(interact[[int]]))
      {
        interact_term <- interact_term * data[,interact[[int]][t]]
      }
      data <- cbind(data, interact_term)
      colnames(data)[ncol(data)] <- paste0(interact[[int]][], collapse = "x")
      X_val_interact <- c(X_val_interact, paste0(interact[[int]][], collapse = "x"))
    }
  }
  # ---------------------------------------- Create interaction terms
  
  # Standardize X_val, X_unval and C to N(0,1) ----------------------
  re_data <- data
  if(rescale)
  {
    ## Shift by the sample mean ---------------------------------------
    re_shift <- c(0, as.numeric(colMeans(data[,c(X_val, X_val_interact, X_unval, C)], na.rm = TRUE)))
    ## Scale inversely by sample standard deviation -------------------
    re_scale <- c(1, as.numeric(apply(data[,c(X_val, X_val_interact, X_unval, C)], MARGIN = 2, FUN = sd, na.rm = TRUE)))
    for (p in 1:length(c(X_val, X_val_interact, X_unval, C)))
    {
      re_data[,c(X_val, X_val_interact, X_unval, C)[p]] <- (re_data[,c(X_val, X_val_interact, X_unval, C)[p]] - re_shift[p+1])/re_scale[p+1]
    }
  }
  # ---------------------- Standardize X_val, X_unval and C to N(0,1)
  
  # Save distinct X -------------------------------------------------
  x_obs <- data.frame(unique(re_data[1:n_v,c(X_val, X_val_interact)]))
  x_obs <- data.frame(x_obs[order(x_obs[,1]),])
  m <- nrow(x_obs)
  x_obs_stacked <- do.call(rbind, replicate(n = (n-n_v), expr = x_obs, simplify = FALSE))
  x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[,1]),])
  colnames(x_obs) <- colnames(x_obs_stacked) <- c(X_val, X_val_interact)
  x_interact_obs <- x_obs[,X_val_interact]
  x_interact_obs_stacked <- x_obs_stacked[,X_val_interact]
  
  # Save static (X*,Y*,X,Y,C) since they don't change ---------------
  comp_dat_val <- re_data[c(1:n_v),c(Y_unval, X_unval, C, Bspline, X_val, X_val_interact, Y_val)]
  comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
  comp_dat_val <- comp_dat_val[,c(Y_unval, X_unval, C, Bspline, X_val, X_val_interact, Y_val, "k")]
  comp_dat_val <- data.matrix(comp_dat_val)
  # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
  # one row per x) --------------------------------------------------
  suppressWarnings(comp_dat_unval <- cbind(re_data[-c(1:n_v),c(Y_unval, X_unval, C, Bspline)],
                                           x_obs_stacked))
  comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
  comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
  colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
  comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1), 
                          k = rep(rep(seq(1,m), each = (n-n_v)), times = 2)))
  
  comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  
  theta_formula <- as.formula(paste0(Y_val, "~", paste(c(X_val, gsub("x", "*", X_val_interact), C), collapse = "+")))
  theta_design_mat <- cbind(int = 1, comp_dat_all[,c(X_val,X_val_interact,C)])

  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[,c(X_unval, Y_val, X_val, C)])
  
  # Standardize Y_val in gamma_design_mat to N(0,1) -----------------
  ## Shift by the sample mean ---------------------------------------
  ## Scale inversely by sample standard deviation -------------------
  ## Create artificially scaled data set ----------------------------
  if(rescale)
  {
    gamma_design_mat[,Y_val] <- (gamma_design_mat[,Y_val] - mean(gamma_design_mat[,Y_val]))/sd(gamma_design_mat[,Y_val])  
  }
  ## Everything else was already standardized above -----------------
  # -------------------------------  Standardize Y_val in gamma_design_mat to N(0,1)
  
  # If unvalidated variable was left blank, assume error-free -------
  ## Need to write simplification here 
  
  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data")))
  {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params <- "Zero"
  }
  if(initial_lr_params == "Zero")
  {
    num_pred <- length(X_val) + length(X_val_interact) + length(C) #preds in analysis model --
    prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
    prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)
  }
  if(initial_lr_params == "Complete-data")
  {
    prev_theta <- theta0 <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(re_data[c(1:n_v),]))$coefficients, ncol = 1)
    prev_gamma <- gamma0 <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(re_data[c(1:n_v),]))$coefficient, ncol = 1)
  }
  # if(initial_lr_params == "Naive")
  # {
  #   prev_theta <- theta0 <- matrix(glm(formula = as.formula(paste0(Y_unval, "~", paste(c(X_unval, C), collapse = "+"))), family = "binomial", data = data.frame(rescaled_data[c(1:n_v),]))$coefficients, ncol = 1)
  #   prev_gamma <- theta0 <- matrix(glm(formula = as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+"))), family = "binomial", data = data.frame(rescaled_data[c(1:n_v),]))$coefficient, ncol = 1)
  # }
  
  # Initialize B-spline coefficients {p_kj}  ------------ 
  ## Numerators sum B(Xi*) over k = 1,...,m -------------
  ## Save as p_val_num for updates ----------------------
  ## (contributions don't change) -----------------------
  p_val_num <- rowsum(x = comp_dat_val[,Bspline], group = comp_dat_val[,"k"], reorder = TRUE)
  prev_p <- p0 <-  t(t(p_val_num)/colSums(p_val_num))
  
  y0 <- comp_dat_unval[,Y_val] == 0
  ystar0 <- comp_dat_unval[,Y_unval] == 0
  
  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1
  
  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric(theta_design_mat[-c(1:n_v),] %*% prev_theta)
    pY_X <- 1/(1 + exp(-mu_theta)) 
    pY_X[y0] <- 1-pY_X[y0]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[-c(1:n_v),] %*% prev_gamma)))
    pYstar[ystar0] <- 1 - pYstar[ystar0]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    pX <- prev_p[rep(rep(seq(1,m), each = (n-n_v)), times = 2),]
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
    #### Avoid NaN resulting from dividing by 0 -----------------------
    psi_t_num_sumover_k_y_j[psi_t_num_sumover_k_y_j == 0] <- 1
    ### And divide them! ----------------------------------------------
    psi_t <- psi_t_num/psi_t_num_sumover_k_y_j
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
    post_multiply <- c((exp(-mu)/(1+exp(-mu)))*(exp(-mu)/(1+exp(-mu))-1))*w_t*theta_design_mat
    hessian_theta <- apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    
    new_theta <- tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta,
                          error = function(err) {matrix(NA, nrow = nrow(prev_theta))
                          })
    if (TRUE %in% is.na(new_theta))
    {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    #if(VERBOSE) print(new_theta)
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
    #if(VERBOSE) print(new_gamma)
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma)<TOL
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., n -----------
    new_p_num <- p_val_num + 
      rowsum(u_t, group = rep(seq(1,m), each = (n-n_v)), reorder = TRUE)
    new_p <- t(t(new_p_num)/colSums(new_p_num))
    #if(VERBOSE) print(new_p[1,])
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p)<TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    
    if(VERBOSE & it%%25 == 0) print(paste("Iteration", it, "complete."))
    
    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    
    it <- it + 1
    
    max_tol_theta = max(abs(new_theta - prev_theta))
    max_tol_gamma = max(abs(new_gamma - prev_gamma))
    max_tol_p = max(abs(new_p - prev_p))
    
    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p 
  }
  
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
                #max_tol_theta = max_tol_theta,
                #max_tol_gamma = max_tol_gamma,
                #max_tol_p = max_tol_p,
                #gamma_at_conv = re_gamma, 
                #p_at_conv = new_p,
                od_loglik_at_conv = NA,
                initial_vals = initial_lr_params))
  }
  if(CONVERGED)
  {
    CONVERGED_MSG <- "Converged"
    #print(paste("SMLE converged after", it, "iterations", sep = " "))
  }
  # ---------------------------------------------- Estimate theta using EM
  if(noSE)
  {
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(n = n, n_v = n_v,
                                            Y_unval=Y_unval, Y_val=Y_val,
                                            X_unval=X_unval, X_val=X_val,
                                            X_val_interact = X_val_interact,
                                            C=C, Bspline=Bspline,
                                            comp_dat_all = comp_dat_all,
                                            theta = new_theta,
                                            gamma = new_gamma,
                                            p = new_p)
    
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
                p_at_conv = new_p,
                od_loglik_at_conv = od_loglik_theta,
                initial_vals = initial_lr_params))
  } else
  {
    rownames(new_theta) <- c("Intercept", colnames(theta_design_mat)[-1])
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_n <- h_n_scale*n^(-1/2) # perturbation ----------------------------
    
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(n = n, n_v = n_v, 
                                            Y_unval=Y_unval, Y_val=Y_val, 
                                            X_unval=X_unval, X_val=X_val, 
                                            X_val_interact = X_val_interact,
                                            C=C, Bspline=Bspline, 
                                            comp_dat_all = comp_dat_all, 
                                            theta = new_theta,
                                            gamma = new_gamma, 
                                            p = new_p)
    
    I_theta <- matrix(od_loglik_theta, nrow = nrow(new_theta), ncol = nrow(new_theta))
    
    #single perts: 271.658 sec elapsed
    single_pert_theta <- sapply(X = seq(1,ncol(I_theta)), FUN = pl_theta, theta = new_theta, h_n = h_n, n_v = n_v, n = n, 
                                       Y_unval=Y_unval, Y_val=Y_val, X_unval=X_unval, X_val=X_val, X_val_interact = X_val_interact, C=C, Bspline=Bspline, 
                                       comp_dat_all = comp_dat_all, gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num)
    
    spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)), ncol = ncol(I_theta), byrow = FALSE)
    #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
    I_theta <- I_theta - spt_wide - t(spt_wide)
    
    for (c in 1:ncol(I_theta))
    {
      pert_theta <- new_theta; pert_theta[c] <- pert_theta[c] + h_n
      ## Avoid 
      double_pert_theta <- sapply(X = seq(c,ncol(I_theta)), FUN = pl_theta, theta = pert_theta, h_n = h_n, n_v = n_v, n = n, 
                                  Y_unval=Y_unval, Y_val=Y_val, X_unval=X_unval, X_val=X_val, X_val_interact = X_val_interact, C=C, Bspline=Bspline, 
                                  comp_dat_all = comp_dat_all, gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num)
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
                p_at_conv = new_p,
                od_loglik_at_conv = od_loglik_theta,
                initial_vals = initial_lr_params))
  }
}