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
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param gamma0 Starting values for `gamma`, the parameters for the outcome error model (a column vector)
#' @param p0 Starting values for `p`, the B-spline coefficients for the approximated covariate error model (a matrix)
#' @param p_val_num Contributions of validated subjects to the numerator for `p`, which are fixed (a matrix)
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return Profile likelihood for `theta`: the value of the observed-data log-likelihood after profiling out other parameters.

profile_out_rw <- function(theta, n, N, Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL, Bspline = NULL,
                           comp_dat_all, theta_pred, gamma_pred, gamma0, p0, p_val_num, TOL, MAX_ITER) {
  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsY <- errorsX <- TRUE
  if (is.null(Y_unval)) {errorsY <- FALSE}
  if (is.null(X_unval) & is.null(X_val)) {errorsX <- FALSE}
  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  sn <- ncol(p0)
  m <- nrow(p0)

  prev_gamma <- gamma0
  prev_p <- p0

  comp_dat_unval <- comp_dat_all[- c(1:n), ]

  theta_design_mat <- cbind(int = 1, comp_dat_all[, theta_pred])

  ### P(Y|X) --------------------------------------------------------
  I_y0 <-  comp_dat_unval[, Y_val] == 0
  mu_theta <- as.numeric(theta_design_mat[- c(1:n), ] %*% theta)
  pY_X <- 1 / (1 + exp(- mu_theta))
  pY_X[I_y0] <- 1 - pY_X[I_y0]
  ### -------------------------------------------------------- P(Y|X)

  if (errorsY) {
    gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
    gamma_design_mat <- cbind(int = 1, comp_dat_all[, gamma_pred])
  }

  CONVERGED <- FALSE
  CONVERGED_MSG <- "SE not converged"
  it <- 1
  # Estimate gamma/p using EM ----------------------------------------------
  while(it < MAX_ITER & !CONVERGED) {
    ###################################################################
    # E Step ----------------------------------------------------------
    ###################################################################
    ## P(Y*|X*,Y,X) ---------------------------------------------------
    if (errorsY) {
      I_ystar0 <- comp_dat_unval[, Y_unval] == 0
      pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[- c(1:n), ] %*% prev_gamma)))
      pYstar[I_ystar0] <- 1 - pYstar[I_ystar0]

    } else { pYstar <- rep(1, times = nrow(comp_dat_unval)) }
    ## --------------------------------------------------- P(Y*|X*,Y,X)
    ## P(X|X*) --------------------------------------------------------
    if (errorsX) {
      if (errorsY) {
        ### p_kj x the B-spline terms
        pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ] * comp_dat_unval[, Bspline]
      } else {
        ### p_kj x the B-spline terms
        pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 1), ] * comp_dat_unval[, Bspline]
      }
    } else { pX <- rep(1, times = nrow(comp_dat_unval)) }
    ## -------------------------------------------------------- P(X|X*)
    ###################################################################
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### Update numerator: P(Y|X,C)*P(Y*|X*,Y,X,C)p_kjB(X*) ------------
    psi_num <- pY_X * pYstar * pX
    ### Update denominator --------------------------------------------
    ### (1) Sum up all rows per id (e.g. sum over xk/y) ---------------
    if (errorsY & errorsX) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2 * m))
    } else if (errorsX) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
    } else if (errorsY) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2))
    }
    ### (2) Then sum over the sn splines ------------------------------
    psi_denom <- rowSums(psi_denom)
    ### Avoid NaN resulting from dividing by 0 -----------------------
    psi_denom[psi_denom == 0] <- 1
    ### And divide them! ----------------------------------------------
    psi_t <- psi_num / psi_denom
    ## ------------------- Update the psi_kyji for unvalidated subjects
    if (errorsX) {
      ## Update the w_kyi for unvalidated subjects --------------------
      ## by summing across the splines/ columns of psi_t --------------
      ## w_t is ordered by i = (N-n), ..., N --------------------------
      w_t <- rowSums(psi_t)
      ## -------------------- Update the w_kyi for unvalidated subjects
      if (errorsY) {
        ## Update the u_kji for unvalidated subjects ------------------
        ## by summing over Y = 0/1 w/i each i, k ----------------------
        ## add top half of psi_t (y = 0) to bottom half (y = 1) -------
        u_t <- psi_t[c(1:(m * (N - n))), ] + psi_t[-c(1:(m * (N - n))), ]
        ## ------------------ Update the u_kji for unvalidated subjects
      } else { u_t <- psi_t }
    } else if (errorsY) {
      w_t <- psi_t
      ## Update the u_kji for unvalidated subjects --------------------
      ## by summing over Y = 0/1 w/i each i, k ------------------------
      ## add top half of psi_t (y = 0) to bottom half (y = 1) ---------
      # u_t <- psi_t[c(1:(N - n))] + psi_t[-c(1:(N - n))]
      ## -------------------- Update the u_kji for unvalidated subjects
    }
    ###################################################################
    # E Step ----------------------------------------------------------
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    if (errorsY) {
      w_t <- c(rep(1, n), w_t)
      mu <- gamma_design_mat %*% prev_gamma
      ### Gradient ----------------------------------------------------
      gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(-mu)/(1 + exp(-mu)))) * gamma_design_mat)), ncol = 1)
      ### ---------------------------------------------------- Gradient
      ### Hessian -----------------------------------------------------
      post_multiply <- c(w_t * (exp(-mu) / (1 + exp(-mu))) * (exp(-mu) / (1 + exp(-mu)) - 1)) * gamma_design_mat
      hessian_gamma <- apply(gamma_design_mat,
                             MARGIN = 2,
                             FUN = hessian_row,
                             pm = post_multiply)
      ### ----------------------------------------------------- Hessian
      ### Newton step -------------------------------------------------
      new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
                            error = function(err) {
                              matrix(NA, nrow = nrow(prev_gamma))
                            })
      if (any(is.na(new_gamma))) {
        suppressWarnings(
          new_gamma <- matrix(glm(formula = gamma_formula,
                                  family = "binomial",
                                  data = data.frame(comp_dat_all),
                                  weights = w_t)$coefficients, ncol = 1)
        )
      }
      ### ------------------------------------------------- Newton step
      ### Check gamma for convergence ---------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
      ### --------------------------------- Check gamma for convergence
    } else { gamma_conv <- TRUE }
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    if (errorsX) {
      ### Update numerators by summing u_t over i = 1, ..., N ---------
      new_p_num <- p_val_num +
        rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
      ### --------- Update numerators by summing u_t over i = 1, ..., N
      ### Standardize so that p_kj sums to 1 over j -------------------
      new_p <- t(t(new_p_num) / colSums(new_p_num))
      ### ------------------- Standardize so that p_kj sums to 1 over j
      ### Check p for convergence -------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
      ### ------------------------------------- Check p for convergence
    } else { p_conv <- TRUE }
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################

    # Check gamma and {p_kj} for overall EM convergence ---------------
    all_conv <- c(gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE
    # --------------- Check gamma and {p_kj} for overall EM convergence

    # Update values for next iteration --------------------------------
    it <- it + 1
    if (errorsY) { prev_gamma <- new_gamma }
    if (errorsX) { prev_p <- new_p }
    # -------------------------------- Update values for next iteration
  }
  if(it == MAX_ITER & !CONVERGED) {
    CONVERGED_MSG <- "MAX_ITER reached"
    new_gamma <- matrix(NA, nrow = nrow(gamma0), ncol = 1)
    new_p <- matrix(NA, nrow = nrow(p0), ncol = ncol(p0))
  }
  if(CONVERGED) CONVERGED_MSG <- "converged"
  if (!errorsY) { new_gamma <- NA }
  if (!errorsX) { new_p <- NA }
  # ---------------------------------------------- Estimate theta using EM
  return(list("psi_at_conv" = psi_t,
              "gamma_at_conv" = new_gamma,
              "p_at_conv" = new_p,
              "converged" = CONVERGED,
              "converged_msg" = CONVERGED_MSG))
}
