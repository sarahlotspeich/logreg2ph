sensY <- 0.95
specY <- 0.90
sim_seed <- sample(seq(1,10000), 1)

prevZ <- 0.25
num_sims <- 101

beta0 <- -1
beta1 <- 0.01
beta2 <- -0.1

beta0_Xb <- -0.65
beta1_Xb <- -0.20
beta2_Xb <- -0.05

theta0 <- -log(specY/(1-specY))
theta1 <- -theta0 - log((1-sensY)/sensY)
theta2 <- 0.01
theta3 <- 0.01
theta4 <- -0.2

N <- 1000
pv <- 0.1
n <- pv*N

Z <- rbinom(n = N, size = 1, prob = prevZ)
X <- rnorm(n = N, mean = 25, sd = 1)
Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(beta0 + beta1*X + beta2*Z)))^(-1))
#Xb <- as.numeric(X >= 25)
#Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(beta0_Xb + beta1_Xb*Xb + beta2_Xb*Z)))^(-1))
U <- rnorm(n = N, mean = 0, sd = 1)
Xstar <- X + U
#Xstarb <- as.numeric(Xstar >= 25)
Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(-(theta0 + theta1*Y + theta2*X + theta3*Xstar + theta4*Z)))^(-1))

V <- seq(1,N) %in% sample(x = seq(1,N), size = n, replace = FALSE)

nsieve <- 20
B <- matrix(0, nrow = N, ncol = nsieve)
B[which(Z == 0),1:15] <- splines::bs(x = Xstar[which(Z == 0)], df = 15, Boundary.knots = range(Xstar[which(Z == 0)]), intercept = TRUE)
B[which(Z == 1),16:20] <- splines::bs(x = Xstar[which(Z == 1)], df = 5, Boundary.knots = range(Xstar[which(Z == 1)]), intercept = TRUE)
colnames(B) <- paste0("bs", seq(1, nsieve))

sdat <- cbind(Y, X, Ystar, Xstar, Z, V, B)
sdat[!V,c("Y", "X")] <- NA

Y_unval = "Ystar"
Y_val = "Y"
X_unval = "Xstar"
X_val = "X"
C = "Z"
Validated = "V"
Bspline = colnames(B)
data = sdat
noSE = FALSE
VERBOSE = FALSE
MAX_ITER = 5000
initial_lr_params = "Zero"
h_n_scale = 1
TOL = 1E-4

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
              h_n = h_n,
              converged = FALSE,
              converged_msg = "Empty sieve in validated data.",
              initial_vals = NA, 
              iterations = 0))
}
# ------------------------------------------ Add the B spline basis

# Standardize X_val and C to N(0,1) -------------------------------
## Shift by the sample mean ---------------------------------------
re_shift <- c(0, as.numeric(colMeans(data[,c(X_val, C)], na.rm = TRUE)))
## Scale inversely by sample standard deviation -------------------
re_scale <- c(1, as.numeric(apply(data[,c(X_val, C)], MARGIN = 2, FUN = sd, na.rm = TRUE)))
## Create artificially scaled data set ----------------------------
re_data <- data
for (p in 1:length(c(X_val, C)))
{
  re_data[,c(X_val, C)[p]] <- (re_data[,c(X_val, C)[p]] - re_shift[p+1])/re_scale[p+1]
}
# ------------------------------- Standardize X_val and C to N(0,1)

# Save distinct X -------------------------------------------------
x_obs <- data.frame(unique(re_data[1:n_v,c(X_val)]))
x_obs <- data.frame(x_obs[order(x_obs[,1]),])
colnames(x_obs) <- c(X_val)
m <- nrow(x_obs)
x_obs_stacked <- do.call(rbind, replicate(n = (n-n_v), expr = x_obs, simplify = FALSE))
x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[,1]),])
colnames(x_obs_stacked) <- c(X_val)
#suppressMessages(data %<>% dplyr::left_join(data.frame(x_k, k = 1:m)))

# Save static (X*,Y*,X,Y,C) since they don't change ---------------
comp_dat_val <- re_data[c(1:n_v),c(Y_unval, X_unval, C, Bspline, X_val, Y_val)]
comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
comp_dat_val <- comp_dat_val[,c(Y_unval, X_unval, C, Bspline, X_val, Y_val, "k")]
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
# Initialize parameter values -------------------------------------
## theta, gamma ---------------------------------------------------
if(!(initial_lr_params %in% c("Zero", "Complete-data", "Naive")))
{
  initial_lr_params <- "Zero"
}
if(initial_lr_params == "Zero")
{
  num_pred <- length(X_val) + length(C) #preds in analysis model --
  prev_theta <- theta0 <- matrix(0, nrow = (num_pred+1), ncol = 1)
  prev_gamma <- gamma0 <- matrix(0, nrow = (length(Y_val) + length(X_unval) + num_pred + 1), ncol = 1)
}
if(initial_lr_params == "Complete-data")
{
  prev_theta <- theta0 <- matrix(glm(formula = as.formula(paste0(Y_val, "~", paste(c(X_val, C), collapse = "+"))), family = "binomial", data = data.frame(rescaled_data[c(1:n_v),]))$coefficients, ncol = 1)
  prev_gamma <- theta0 <- matrix(glm(formula = as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+"))), family = "binomial", data = data.frame(rescaled_data[c(1:n_v),]))$coefficient, ncol = 1)
}
if(initial_lr_params == "Naive")
{
  prev_theta <- theta0 <- matrix(glm(formula = as.formula(paste0(Y_unval, "~", paste(c(X_unval, C), collapse = "+"))), family = "binomial", data = data.frame(rescaled_data[c(1:n_v),]))$coefficients, ncol = 1)
  prev_gamma <- theta0 <- matrix(glm(formula = as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+"))), family = "binomial", data = data.frame(rescaled_data[c(1:n_v),]))$coefficient, ncol = 1)
}

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
  #### Avoid NaN resulting from dividing by 0 -----------------------
  psi_t_num_sumover_k_y_j[psi_t_num_sumover_k_y_j == 0] <- 1
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
  
  new_theta <- tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta,
                        error = function(err) {matrix(NA, nrow = nrow(prev_theta))
                        })
  if (TRUE %in% is.na(new_theta))
  {
    suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
  }
  if(VERBOSE) print(new_theta)
  ### Check for convergence -----------------------------------------
  theta_conv <- abs(new_theta - prev_theta)<TOL
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
  
  # Update values for next iteration  -------------------------------
  prev_theta <- new_theta
  prev_gamma <- new_gamma
  prev_p <- new_p 
}

rownames(new_theta) <- c("Intercept", X_val,C)
# Estimate Cov(theta) using profile likelihood -------------------------
h_n <- h_n_scale*n^(-1/2) # perturbation ----------------------------

## Calculate pl(theta) -------------------------------------------------
od_loglik_theta <- observed_data_loglik(n = n, n_v = n_v, 
                                        Y_unval=Y_unval, Y_val=Y_val, 
                                        X_unval=X_unval, X_val=X_val, 
                                        C=C, Bspline=Bspline, 
                                        comp_dat_all = comp_dat_all, 
                                        theta = new_theta,
                                        gamma = new_gamma, 
                                        p = new_p)

I_theta <- matrix(0, nrow = nrow(new_theta), ncol = nrow(new_theta))
for (k in 1:ncol(I_theta))
{
  pert_k <- new_theta; pert_k[k] <- pert_k[k] + h_n
  pl_params <- profile_out(theta = pert_k, 
                           n_v = n_v, n = n, 
                           Y_unval=Y_unval, Y_val=Y_val, 
                           X_unval=X_unval, X_val=X_val, 
                           C=C, Bspline=Bspline, 
                           comp_dat_all = comp_dat_all, 
                           gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num)
  od_loglik_pert_k <- observed_data_loglik(n = n, n_v = n_v, 
                                           Y_unval=Y_unval, Y_val=Y_val, 
                                           X_unval=X_unval, X_val=X_val, 
                                           C=C, Bspline=Bspline, 
                                           comp_dat_all = comp_dat_all, 
                                           theta = pert_k,
                                           gamma = pl_params$gamma, 
                                           p = pl_params$p_at_conv)
  for (l in k:nrow(I_theta))
  {
    pert_l <- new_theta; pert_l[l] <- pert_l[l] + h_n
    pert_both <- pert_l; pert_both[k] <- pert_both[k] + h_n
    
    pl_params <- profile_out(theta = pert_both, 
                             n_v = n_v, n = n, 
                             Y_unval=Y_unval, Y_val=Y_val, 
                             X_unval=X_unval, X_val=X_val, 
                             C=C, Bspline=Bspline, 
                             comp_dat_all = comp_dat_all, 
                             gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num)
    
    od_loglik_pert_both <- observed_data_loglik(n = n, n_v = n_v, 
                                                Y_unval=Y_unval, Y_val=Y_val, 
                                                X_unval=X_unval, X_val=X_val, 
                                                C=C, Bspline=Bspline, 
                                                comp_dat_all = comp_dat_all, 
                                                theta = pert_both, 
                                                gamma = pl_params$gamma, 
                                                p = pl_params$p_at_conv)
    
    if (l == k)
    {
      I_theta[l,k] <- od_loglik_pert_both - 2*od_loglik_pert_k + od_loglik_theta
    } else
    {
      pl_params <- profile_out(theta = pert_l, 
                               n_v = n_v, n = n, 
                               Y_unval=Y_unval, Y_val=Y_val, 
                               X_unval=X_unval, X_val=X_val, 
                               C=C, Bspline=Bspline, 
                               comp_dat_all = comp_dat_all, 
                               gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num)
      
      od_loglik_pert_l <- observed_data_loglik(n = n, n_v = n_v, 
                                               Y_unval=Y_unval, Y_val=Y_val, 
                                               X_unval=X_unval, X_val=X_val, 
                                               C=C, Bspline=Bspline, 
                                               comp_dat_all = comp_dat_all, 
                                               theta = pert_l, 
                                               gamma = pl_params$gamma, 
                                               p =  pl_params$p_at_conv)
      I_theta[k,l] <- I_theta[l,k] <- od_loglik_pert_both - od_loglik_pert_k - od_loglik_pert_l + od_loglik_theta # symmetry of covariance 
    }
  }
}
I_theta <- h_n^(-2) * I_theta
cov_theta <- -solve(I_theta)
# ------------------------- Estimate Cov(theta) using profile likelihood
re_theta <- new_theta
re_theta[-1] <- re_theta[-1]/re_scale[-1]
re_theta[1] <- re_theta[1] - sum(re_theta[-1]*re_shift[-1])

re_se_theta <- sqrt(diag(cov_theta))
re_se_theta[-1] <- re_se_theta[-1]/re_scale[-1]
re_se_theta[1] <- cov_theta[1,1] + sum(diag(cov_theta)[-1]*(re_shift[-1]/re_scale[-1])^2)
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

list(Coefficients = data.frame(Coefficient = re_theta, 
                                      SE = re_se_theta),
            h_n = h_n,
            converged = CONVERGED,
            converged_msg = CONVERGED_MSG,
            initial_vals = initial_lr_params, 
            iterations = it)
