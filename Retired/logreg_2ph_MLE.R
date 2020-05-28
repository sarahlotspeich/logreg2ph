setup_complete_data <- function(dat, Y_unval, Y_val, contX_unval, contX_val, 
                                binX_unval, binX_val, C = NULL, Validated)
{
  n <- nrow(dat)
  n_v <- sum(dat[,Validated])
  
  # Reorder so that the n_v validated subjects are first ------------
  dat <- dat[order(as.numeric(dat[,Validated]), decreasing = TRUE),]
  
  # Save distinct X -------------------------------------------------
  x_obs <- data.frame(unique(dat[1:n_v,c(contX_val, binX_val)]))
  x_obs <- data.frame(x_obs[order(x_obs[,1]),])
  colnames(x_obs) <- c(contX_val, binX_val)
  m <- nrow(x_obs)
  x_obs_stacked <- do.call(rbind, replicate(n = (n-n_v), expr = x_obs, simplify = FALSE))
  x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[,1]),])
  colnames(x_obs_stacked) <- c(contX_val, binX_val)
  
  # Save static (X*,Y*,X,Y,C) since they don't change ---------------
  comp_dat_val <- dat[c(1:n_v),c(Y_unval, contX_unval, binX_unval, C, contX_val, binX_val, Y_val)]
  comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
  comp_dat_val <- comp_dat_val[,c(Y_unval, contX_unval, binX_unval, C, contX_val, binX_val, Y_val, "k")]
  comp_dat_val <- data.matrix(comp_dat_val)
  # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
  # one row per x) --------------------------------------------------
  suppressWarnings(comp_dat_unval <- cbind(dat[-c(1:n_v),c(Y_unval, contX_unval, binX_unval, C)],
                                           x_obs_stacked))
  comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
  comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
  colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
  comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1), 
                                      k = rep(rep(seq(1,m), each = (n-n_v)), times = 2)))
  
  comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  return(comp_dat_all)
}

observed_data_loglik_fp <- function(params, n, n_v, m, Y_unval=NULL, Y_val=NULL, 
                                 contX_unval=NULL, contX_val=NULL, binX_unval=NULL, binX_val=NULL, C=NULL, 
                                 comp_dat_all)
{
  X_val <- c(contX_val, binX_val)
  X_unval <- c(contX_unval, binX_unval)
  
  theta <- params[1:(length(X_val)+length(C)+1)]
  d_theta <- length(theta)
  gamma <- params[c(1:(length(X_val)+length(X_unval)+length(C)+2))+d_theta]
  d_gamma <- length(gamma)
  sigma <- params[length(params)]
  lambda <- params[-c(1:(d_theta+d_gamma), length(params))]
  d_lambda <- length(lambda)
  
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
  ## Sum over log[P(Xi|Xi*)] ------------------------------------------------------
  mu <- cbind(int = 1, comp_dat_all[c(1:n_v),c(binX_val, C, contX_unval, binX_unval)]) %*% lambda[1:6]
  pX <- 1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(comp_dat_all[c(1:n_v), contX_val] - mu)^2)
  return_loglik <- return_loglik + sum(log(pX))
  
  mu <- cbind(int = 1, comp_dat_all[c(1:n_v),c(binX_val[2], C, contX_unval, binX_unval)]) %*% lambda[7:11]
  pX <- 1/(1 + exp(-mu))
  pX <- ifelse(as.vector(comp_dat_all[c(1:n_v),binX_val[1]]) == 0, 1-pX, pX)
  return_loglik <- return_loglik + sum(log(pX))
  
  mu <- cbind(int = 1, comp_dat_all[c(1:n_v),c(C, contX_unval, binX_unval)]) %*% lambda[12:15]
  pX <- 1/(1 + exp(-mu))
  pX <- ifelse(as.vector(comp_dat_all[c(1:n_v),binX_val[2]]) == 0, 1-pX, pX)
  return_loglik <- return_loglik + sum(log(pX))
  ## ------------------------------------------------------ Sum over log[P(Xi|Xi*)]
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
  ## Sum over log[P(Xi|Xi*)] ------------------------------------------------------
  mu <- cbind(int = 1, comp_dat_all[-c(1:n_v),c(binX_val, C, contX_unval, binX_unval)]) %*% lambda[1:6]
  pX1 <- 1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(comp_dat_all[-c(1:n_v), contX_val] - mu)^2)
  
  mu <- cbind(int = 1, comp_dat_all[-c(1:n_v),c(binX_val[2], C, contX_unval, binX_unval)]) %*% lambda[7:11]
  pX2 <- 1/(1 + exp(-mu))
  pX2 <- ifelse(as.vector(comp_dat_all[-c(1:n_v),binX_val[1]]) == 0, 1-pX2, pX2)
  
  mu <- cbind(int = 1, comp_dat_all[-c(1:n_v),c(C, contX_unval, binX_unval)]) %*% lambda[12:15]
  pX3 <- 1/(1 + exp(-mu))
  pX3 <- ifelse(as.vector(comp_dat_all[-c(1:n_v),binX_val[2]]) == 0, 1-pX3, pX3)
  ## ------------------------------------------------------ Sum over log[P(Xi|Xi*)]
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x p_kj ---------------------------
  person_sum <- rowsum(pY_X*pYstar*pX1*pX2*pX3, group = rep(seq(1,(n-n_v)), times = 2*m))
  log_person_sum <- log(person_sum)
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(-return_loglik)
}


gradient_logreg_2ph_MLE <- function(theta, gamma, lambda1, lambda2, lambda3, sigma, 
                                    n, m, complete_data, 
                                    Y_unval, Y_val, contX_unval, contX_val, binX_unval, binX_val)
{
  mu1_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_val, binX_val)])) %*% theta
  mu2_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)])) %*% gamma
  mu3_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(binX_val, contX_unval, binX_unval)])) %*% lambda1
  mu4_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(binX_val[2], contX_unval, binX_unval)])) %*% lambda2
  mu5_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_unval, binX_unval)])) %*% lambda3
  
  pY_unval <- exp(-(1-complete_data[-c(1:n),Y_val])*mu1_unval)/(1+exp(-mu1_unval))
  pYstar_unval <- exp(-(1-complete_data[-c(1:n),Y_unval])*mu2_unval)/(1+exp(-mu2_unval))
  pX1_unval <- 1/sqrt(2*pi*sigma^2)*exp(-1/(2*sigma^2)*(complete_data[-c(1:n),contX_val]-mu3_unval)^2)
  pX2_unval <- exp(-(1-complete_data[-c(1:n),binX_val[1]])*mu4_unval)/(1+exp(-mu4_unval))
  pX3_unval <- exp(-(1-complete_data[-c(1:n),binX_val[2]])*mu5_unval)/(1+exp(-mu5_unval))
  joint_unval <- pY_unval * pYstar_unval * pX1_unval * pX2_unval * pX3_unval
  
  delta <- rowsum(x = joint_unval, group = rep(seq(1,(nrow(sim_data)-n)), times = 2*m), reorder = FALSE)
  delta <- rep(delta, times = 2*m)
  
  #### P(Y|X) ---------------------------------------------------------------------------------------------
  mu1_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_val, binX_val)])) %*% theta
  sig1_val <- 1/(1 + exp(-mu1_val))
  pY_val <- (complete_data[c(1:n),Y_val])*sig1_val + (1-complete_data[c(1:n),Y_val])*(1-sig1_val)
  x1_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_val, binX_val)]))
  d_sig1_val <- x1_val * c(sig1_val*(1-sig1_val))
  d_pY_val <- (2*complete_data[c(1:n),Y_val]-1) * d_sig1_val
  grad_ll_theta <- colSums(
    d_pY_val/matrix(data = pY_val, nrow = length(pY_val), ncol = ncol(d_pY_val), byrow = FALSE)
  )
  
  x1_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_val, binX_val)]))
  mu1_unval <- x1_val %*% theta
  sig1_unval <- 1/(1 + exp(-mu1_unval))
  d_sig1_unval <- x1_unval * c(sig1_unval*(1-sig1_unval))
  d_pY_unval <- (2*complete_data[-c(1:n),Y_val]-1) * d_sig1_unval
  grad_ll_theta <- grad_ll_theta + 
    colSums(
      d_pY_unval * c(pYstar_unval * pX1_unval * pX2_unval * pX3_unval / delta)
    )
  #### --------------------------------------------------------------------------------------------- P(Y|X)
  
  #### P(Y*|X*,Y,X) ---------------------------------------------------------------------------------------
  mu2_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)])) %*% gamma
  sig2_val <- 1/(1 + exp(-mu2_val))
  pYstar_val <- (complete_data[c(1:n),Y_unval])*sig2_val + (1-complete_data[c(1:n),Y_unval])*(1-sig2_val)
  x2_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)]))
  d_sig2_val <- x2_val * c(sig2_val*(1-sig2_val))
  d_pYstar_val <- (2*complete_data[c(1:n),Y_unval]-1) * d_sig2_val
  grad_ll_gamma <- colSums(
    d_pYstar_val/matrix(data = pYstar_val, nrow = length(pYstar_val), ncol = ncol(d_pYstar_val), byrow = FALSE)
  )
  
  x2_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)]))
  sig2_unval <- 1/(1 + exp(-mu2_unval))
  d_sig2_unval <- x2_unval * c(sig2_unval*(1-sig2_unval))
  d_pYstar_unval <- (2*complete_data[-c(1:n),Y_unval]-1) * d_sig2_unval
  grad_ll_gamma <- grad_ll_gamma + 
    colSums(
      d_pYstar_unval * c(pY_unval * pX1_unval * pX2_unval * pX3_unval / delta)
    )
  #### --------------------------------------------------------------------------------------- P(Y*|X*,Y,X)
  
  #### P(X1|X2,X3,X*) -------------------------------------------------------------------------------------
  ##### Vector of parameters: lambda1 ---------------------------------------------------------------------
  x3_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(binX_val,contX_unval, binX_unval)]))
  mu3_val <- x3_val %*% lambda1
  grad_ll_lambda1 <- colSums(x3_val * c((1/sigma^2) * (complete_data[c(1:n),contX_val] - mu3_val)))
  
  x3_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(binX_val,contX_unval, binX_unval)]))
  mu3_unval <- x3_unval %*% lambda1
  d_delta3_unval <- (1/sigma^2)*x3_unval*c((complete_data[-c(1:n),contX_val]-mu3_unval)*joint_unval)
  grad_ll_lambda1 <- grad_ll_lambda1 + colSums(d_delta3_unval/delta)
  ##### Nuisance parameter: sigma -------------------------------------------------------------------------
  grad_ll_sigma <- sum(-(sigma^2-(complete_data[c(1:n), contX_val]-mu3_val)^2)/sigma^3)
  grad_ll_sigma <- grad_ll_sigma + 
    sum(-(sigma^2-(complete_data[-c(1:n), contX_val]-mu3_unval)^2)/sigma^3*joint_unval/delta)
  #### ------------------------------------------------------------------------------------- P(X1|X2,X3,X*)
  
  #### P(X2|X3, X*) ---------------------------------------------------------------------------------------
  mu4_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(binX_val[2], contX_unval, binX_unval)])) %*% lambda2
  sig4_val <- 1/(1 + exp(-mu4_val))
  pX2_val <- (complete_data[c(1:n),binX_val[1]])*sig4_val + (1-complete_data[c(1:n),binX_val[1]])*(1-sig4_val)
  x4_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(binX_val[2], contX_unval, binX_unval)]))
  d_sig4_val <- x4_val * c(sig4_val*(1-sig4_val))
  d_pX2_val <- (2*complete_data[c(1:n),binX_val[1]]-1) * d_sig4_val
  grad_ll_lambda2 <- colSums(
    d_pX2_val/matrix(data = pX2_val, nrow = length(pX2_val), ncol = ncol(d_pX2_val), byrow = FALSE)
  )
  
  x4_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(binX_val[2], contX_unval, binX_unval)]))
  sig4_unval <- 1/(1 + exp(-mu4_unval))
  d_sig4_unval <- x4_unval * c(sig4_unval*(1-sig4_unval))
  d_pX2_unval <- (2*complete_data[-c(1:n),binX_val[1]]-1) * d_sig4_unval
  grad_ll_lambda2 <- grad_ll_lambda2 + 
    colSums(
      d_pX2_unval * c(pY_unval * pYstar_unval * pX1_unval * pX3_unval / delta)
    )
  #### ---------------------------------------------------------------------------------------- P(X2|X3, X*)
  
  #### P(X3|X*) --------------------------------------------------------------------------------------------
  mu5_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval)])) %*% lambda3
  sig5_val <- 1/(1 + exp(-mu5_val))
  pX3_val <- (complete_data[c(1:n),binX_val[2]])*sig5_val + (1-complete_data[c(1:n),binX_val[2]])*(1-sig5_val)
  x5_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval)]))
  d_sig5_val <- x5_val * c(sig5_val*(1-sig5_val))
  d_pX3_val <- (2*complete_data[c(1:n),binX_val[2]]-1) * d_sig5_val
  grad_ll_lambda3 <- colSums(
    d_pX3_val/matrix(data = pX3_val, nrow = length(pX3_val), ncol = ncol(d_pX3_val), byrow = FALSE)
  )
  
  x5_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_unval, binX_unval)]))
  sig5_unval <- 1/(1 + exp(-mu5_unval))
  d_sig5_unval <- x5_unval * c(sig5_unval*(1-sig5_unval))
  d_pX3_unval <- (2*complete_data[-c(1:n),binX_val[2]]-1) * d_sig5_unval
  grad_ll_lambda3 <- grad_ll_lambda3 + 
    colSums(
      d_pX3_unval * c(pY_unval * pYstar_unval * pX1_unval * pX2_unval / delta)
    )
  #### -------------------------------------------------------------------------------------------- P(X3|X*)
  return(matrix(c(grad_ll_theta, grad_ll_gamma, grad_ll_lambda1, grad_ll_lambda2, grad_ll_lambda3, sigma)))
}

hessian_logreg_2ph_MLE <- function(theta, gamma, lambda1, lambda2, lambda3, sigma, 
                                   n, m, complete_data, 
                                   Y_unval, Y_val, contX_unval, contX_val, binX_unval, binX_val)
{
  params <- c(theta, gamma, lambda1, lambda2, lambda3, sigma)
  
  theta_params <- c("Int", contX_val, binX_val)
  gamma_params <- c("Int", contX_unval, binX_unval, Y_val, contX_val, binX_val)
  lambda1_params <- c("Int", binX_val, contX_unval, binX_unval)
  lambda2_params <- c("Int", binX_val[2], contX_unval, binX_unval)
  lambda3_params <- c("Int", contX_unval, binX_unval)
  
  params <- c(theta_params, gamma_params, lambda1_params, lambda2_params, lambda3_params, "Sigma")
  d_H <- length(params)
  H <- H_v <- H_uv <- matrix(data = 0, nrow = d_H, ncol = d_H)
  
  #### P(Y|X) ---------------------------------------------------------------------------------------------
  mu1_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_val, binX_val)])) %*% theta
  sig1_val <- 1/(1 + exp(-mu1_val))
  pY_val <- (complete_data[c(1:n),Y_val])*sig1_val + (1-complete_data[c(1:n),Y_val])*(1-sig1_val)
  x1_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_val, binX_val)]))
  d_sig1_val <- x1_val * c(sig1_val*(1-sig1_val))
  d_pY_val <- (2*complete_data[c(1:n),Y_val]-1) * d_sig1_val
  
  for (p in 1:length(theta_params))
  {
    post_multiply <- c(x1_val[,p]*(2*complete_data[c(1:n),Y_val]-1)*(1-2*sig1_val)*pY_val^(-1)) - d_pY_val[p]*c(pY_val^(-2))
    H_v[p,c(p:length(theta_params))] <- H_v[c(p:length(theta_params)),p] <- colSums(d_sig1_val * post_multiply)[c(p:length(theta_params))]
  }
  
  d2_pY_val <- c((2*complete_data[c(1:n),Y_val]-1) * (1-2*sig1_val) * sig1_val * (1-sig1_val)) * (x1_val)^2
  
  H_theta_val <- matrix(data = 0, nrow = (length(theta_params)+1), ncol = (length(theta_params)+1))
  
  for (r in 1:nrow(H_theta_val))
  {
    for (c in r:ncol(H_theta_val))
    {
      d2_pY_val <- c((2*complete_data[c(1:n),Y_val]-1) * (1-2*sig1_val) * sig1_val * (1-sig1_val) * x1_val[,r]) * x1_val
      H_theta_val[r,c(r:ncol(H_theta_val))] <- H_theta_val[c(r:ncol(H_theta_val)),r] <- colSums(matrix(d2_pY_val[,c(r:ncol(H_theta_val))]) * c(pY_val^(-1)) - matrix(d_pY_val[,c(r:ncol(H_theta_val))]) * d_pY_val[,r] * c(pY_val^(-2)))
    }
  }
  
  
  diag(H_theta_val) <-  colSums(d2_pY_val * c(pY_val^(-1)) - d_pY_val^2 * c(pY_val^(-2)))
  
  grad_ll_theta <- colSums(
    d_pY_val/matrix(data = pY_val, nrow = length(pY_val), ncol = ncol(d_pY_val), byrow = FALSE)
  )
  
  x1_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_val, binX_val)]))
  mu1_unval <- x1_val %*% theta
  sig1_unval <- 1/(1 + exp(-mu1_unval))
  d_sig1_unval <- x1_unval * c(sig1_unval*(1-sig1_unval))
  d_pY_unval <- (2*complete_data[-c(1:n),Y_val]-1) * d_sig1_unval
  d_delta1_unval <- d_pY_unval * c(pYstar_unval * pX1_unval * pX2_unval * pX3_unval)
  grad_ll_theta <- grad_ll_theta + colSums(d_delta1_unval / delta)
  #### --------------------------------------------------------------------------------------------- P(Y|X)
  
  #### P(Y*|X*,Y,X) ---------------------------------------------------------------------------------------
  mu2_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)])) %*% gamma
  sig2_val <- 1/(1 + exp(-mu2_val))
  pYstar_val <- (complete_data[c(1:n),Y_unval])*sig2_val + (1-complete_data[c(1:n),Y_unval])*(1-sig2_val)
  x2_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)]))
  d_sig2_val <- x2_val * c(sig2_val*(1-sig2_val))
  d_pYstar_val <- (2*complete_data[c(1:n),Y_unval]-1) * d_sig2_val
  grad_ll_gamma <- colSums(
    d_pYstar_val/matrix(data = pYstar_val, nrow = length(pYstar_val), ncol = ncol(d_pYstar_val), byrow = FALSE)
  )
  
  x2_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_unval, binX_unval, Y_val, contX_val, binX_val)]))
  sig2_unval <- 1/(1 + exp(-mu2_unval))
  d_sig2_unval <- x2_unval * c(sig2_unval*(1-sig2_unval))
  d_pYstar_unval <- (2*complete_data[-c(1:n),Y_unval]-1) * d_sig2_unval
  d_delta2_unval <- d_pYstar_unval * c(pY_unval * pX1_unval * pX2_unval * pX3_unval)
  grad_ll_gamma <- grad_ll_gamma + colSums(d_delta2_unval / delta)
  #### --------------------------------------------------------------------------------------- P(Y*|X*,Y,X)
  
  #### P(X1|X2,X3,X*) -------------------------------------------------------------------------------------
  ##### Vector of parameters: lambda1 ---------------------------------------------------------------------
  x3_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(binX_val,contX_unval, binX_unval)]))
  mu3_val <- x3_val %*% lambda1
  grad_ll_lambda1 <- colSums(x3_val * c((1/sigma^2) * (complete_data[c(1:n),contX_val] - mu3_val)))
  
  x3_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(binX_val,contX_unval, binX_unval)]))
  mu3_unval <- x3_unval %*% lambda1
  d_delta3_unval <- (1/sigma^2)*x3_unval*c((complete_data[-c(1:n),contX_val]-mu3_unval)*joint_unval)
  grad_ll_lambda1 <- grad_ll_lambda1 + colSums(d_delta3_unval/delta)
  ##### Nuisance parameter: sigma -------------------------------------------------------------------------
  grad_ll_sigma <- sum(-(sigma^2-(complete_data[c(1:n), contX_val]-mu3_val)^2)/sigma^3)
  grad_ll_sigma <- grad_ll_sigma + 
    sum(-(sigma^2-(complete_data[-c(1:n), contX_val]-mu3_unval)^2)/sigma^3*joint_unval/delta)
  #### ------------------------------------------------------------------------------------- P(X1|X2,X3,X*)
  
  #### P(X2|X3, X*) ---------------------------------------------------------------------------------------
  mu4_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(binX_val[2], contX_unval, binX_unval)])) %*% lambda2
  sig4_val <- 1/(1 + exp(-mu4_val))
  pX2_val <- (complete_data[c(1:n),binX_val[1]])*sig4_val + (1-complete_data[c(1:n),binX_val[1]])*(1-sig4_val)
  x4_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(binX_val[2], contX_unval, binX_unval)]))
  d_sig4_val <- x4_val * c(sig4_val*(1-sig4_val))
  d_pX2_val <- (2*complete_data[c(1:n),binX_val[1]]-1) * d_sig4_val
  grad_ll_lambda2 <- colSums(
    d_pX2_val/matrix(data = pX2_val, nrow = length(pX2_val), ncol = ncol(d_pX2_val), byrow = FALSE)
  )
  
  x4_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(binX_val[2], contX_unval, binX_unval)]))
  sig4_unval <- 1/(1 + exp(-mu4_unval))
  d_sig4_unval <- x4_unval * c(sig4_unval*(1-sig4_unval))
  d_pX2_unval <- (2*complete_data[-c(1:n),binX_val[1]]-1) * d_sig4_unval
  d_delta4_unval <- d_pX2_unval * c(pY_unval * pYstar_unval * pX1_unval * pX3_unval)
  grad_ll_lambda2 <- grad_ll_lambda2 + colSums(d_delta4_unval / delta)
  #### ---------------------------------------------------------------------------------------- P(X2|X3, X*)
  
  #### P(X3|X*) --------------------------------------------------------------------------------------------
  mu5_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval)])) %*% lambda3
  sig5_val <- 1/(1 + exp(-mu5_val))
  pX3_val <- (complete_data[c(1:n),binX_val[2]])*sig5_val + (1-complete_data[c(1:n),binX_val[2]])*(1-sig5_val)
  x5_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_unval, binX_unval)]))
  d_sig5_val <- x5_val * c(sig5_val*(1-sig5_val))
  d_pX3_val <- (2*complete_data[c(1:n),binX_val[2]]-1) * d_sig5_val
  grad_ll_lambda3 <- colSums(
    d_pX3_val/matrix(data = pX3_val, nrow = length(pX3_val), ncol = ncol(d_pX3_val), byrow = FALSE)
  )
  
  x5_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_unval, binX_unval)]))
  sig5_unval <- 1/(1 + exp(-mu5_unval))
  d_sig5_unval <- x5_unval * c(sig5_unval*(1-sig5_unval))
  d_pX3_unval <- (2*complete_data[-c(1:n),binX_val[2]]-1) * d_sig5_unval
  d_delta5_unval <- d_pX3_unval * c(pY_unval * pYstar_unval * pX1_unval * pX2_unval)
  grad_ll_lambda3 <- grad_ll_lambda3 + colSums(d_delta5_unval/delta)
  #### -------------------------------------------------------------------------------------------- P(X3|X*)
  
  gradient <- matrix(c(grad_ll_theta, grad_ll_gamma, grad_ll_lambda1, grad_ll_lambda2, grad_ll_lambda3, grad_ll_sigma))
  
  #### P(Y|X) ---------------------------------------------------------------------------------------------
  mu1_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_val, binX_val)])) %*% theta
  sig1_val <- 1/(1 + exp(-mu1_val))
  pY_val <- (complete_data[c(1:n),Y_val])*sig1_val + (1-complete_data[c(1:n),Y_val])*(1-sig1_val)
  x1_val <- data.matrix(cbind(Int = 1, complete_data[c(1:n),c(contX_val, binX_val)]))
  d_sig1_val <- x1_val * c(sig1_val*(1-sig1_val))
  d_pY_val <- (2*complete_data[c(1:n),Y_val]-1) * d_sig1_val
  grad_ll_theta <- colSums(
    d_pY_val/matrix(data = pY_val, nrow = length(pY_val), ncol = ncol(d_pY_val), byrow = FALSE)
  )
  
  x1_unval <- data.matrix(cbind(Int = 1, complete_data[-c(1:n),c(contX_val, binX_val)]))
  mu1_unval <- x1_val %*% theta
  sig1_unval <- 1/(1 + exp(-mu1_unval))
  d_sig1_unval <- x1_unval * c(sig1_unval*(1-sig1_unval))
  d_pY_unval <- (2*complete_data[-c(1:n),Y_val]-1) * d_sig1_unval
  d_delta1_unval <- d_pY_unval * c(pYstar_unval * pX1_unval * pX2_unval * pX3_unval)
  grad_ll_theta <- grad_ll_theta + colSums(d_delta1_unval / delta)
  #### --------------------------------------------------------------------------------------------- P(Y|X)
  
  for (i in 1:d_H)
  {
    
  }
  return(H)
}