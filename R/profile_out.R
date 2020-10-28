#' Profiles out nuisance parameters from the observed-data log-likelihood for a given value of theta
#'
#' For a given vector `theta` to parameterize P(Y|X,C), this function repeats the EM algorithm to find
#' the values of `gamma` and `p` at convergence. The resulting parameters are used to find the profile
#' log-likelihood for `theta` by plugging them into the observed-data log-likelihood.
#' This function is used by `pl_theta()`.
#
#' @param theta Parameters for the analysis model (a column vector)
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index)
#' @param Y_val Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param C (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param gamma0 Starting values for `gamma`, the parameters for the outcome error model (a column vector)
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return Profile likelihood for `theta` after perturbing element `k` by `h_n`.

profile_out <- function(theta, n, N, Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL, Bspline = NULL, comp_dat_all, gamma0, p0, p_val_num, TOL, MAX_ITER) {
  sn <- ncol(p0)
  m <- nrow(p0)

  prev_gamma <- gamma0
  prev_p <- p0

  comp_dat_unval <- comp_dat_all[- c(1:n), ]

  theta_design_mat <- cbind(int = 1, comp_dat_all[, c(X_val, C)])

  ### P(Y|X) --------------------------------------------------------
  y0 <-  comp_dat_unval[, Y_val] == 0
  mu_theta <- as.numeric(theta_design_mat[- c(1:n), ] %*% theta)
  pY_X <- 1 / (1 + exp(- mu_theta))
  pY_X[y0] <- 1 - pY_X[y0]
  ### -------------------------------------------------------- P(Y|X)

  ystar0 <- comp_dat_unval[, Y_unval] == 0
  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[, c(X_unval, Y_val, X_val, C)])

  CONVERGED <- FALSE
  CONVERGED_MSG <- "SE not converged"
  it <- 1
  # Estimate gamma/p using EM ----------------------------------------------
  while(it < MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[- c(1:n), ] %*% prev_gamma)))
    pYstar[ystar0] <- 1 - pYstar[ystar0]
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### p_kj ----------------------------------------------------------
    ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
    pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ]
    ### ---------------------------------------------------------- p_kj
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)p_kjB(X*) ------------------------------
    psi_t_num <- pY_X * pYstar * comp_dat_unval[, Bspline] * pX
    ### Update denominator --------------------------------------------
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = rep(seq(1, (N - n)), times = 2 * m))
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
    u_t <- psi_t[c(1:(m * (N - n))), ] + psi_t[-c(1:(m * (N - n))), ]
    ## ---------------------- Update the u_kji for unvalidated subjects
    ###################################################################
    # E Step ----------------------------------------------------------

    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    w_t <- c(rep(1, n), w_t)
    mu <- gamma_design_mat %*% prev_gamma
    gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(-mu)/(1 + exp(-mu)))) * gamma_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c(w_t * (exp(-mu) / (1 + exp(-mu))) * (exp(-mu) / (1 + exp(-mu)) - 1)) * gamma_design_mat
    hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
                          error = function(err) {
                            matrix(NA, nrow = nrow(prev_gamma))
                          })
    if (any(is.na(new_gamma))) {
      suppressWarnings(
        new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1)
        )
    }
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma) < TOL
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N -----------
    new_p_num <- p_val_num + rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p <- t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p) < TOL
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
  if(it == MAX_ITER & !CONVERGED) {
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
