#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems
#'
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2020)
#
#'
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index). If \code{Y_unval} is null, the outcome is assumed to be error-free.
#' @param Y_val Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index).  If \code{X_unval} and \code{X_val} are \code{null}, all precictors assumed to be error-free.
#' @param X_val Column(s) with the validated predictors (can be name or numeric index). If \code{X_unval} and \code{X_val} are \code{null}, all precictors assumed to be error-free.
#' @param C (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Validated Columns with the validation indicator (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param data A dataframe with one row per subject containing columns: Y_unval, Y_val, X_unval, X_val, C, Validated, and Bspline.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param initial_lr_params Initial values for parametric model parameters. Choices include (1) "Zero" (non-informative starting values) or (2) "Complete-data" (estimated based on validated subjects only)
#' @param h_N_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is `h_N_scale = 1`.
#' @param noSE Indicator for whether standard errors are desired. Defaults to noSE = FALSE.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{Coefficients}{dataframe with final coefficient and standard error estimates (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates.}
#' @export

logreg2ph_rw <- function(Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL,
                      Validated = NULL, Bspline = NULL, data, theta_pred = NULL, gamma_pred = NULL,
                      initial_lr_params = "Zero", h_N_scale = 1, noSE = FALSE, TOL = 1E-4, MAX_ITER = 1000)
{
  N <- nrow(data)
  n <- sum(data[, Validated])

  # Reorder so that the n validated subjects are first ------------
  data <- data[order(as.numeric(data[, Validated]), decreasing = TRUE), ]

  # Determine error setting -----------------------------------------
  ## If unvalidated variable was left blank, assume error-free ------
  errorsY <- errorsX <- TRUE
  if (is.null(Y_unval)) {errorsY <- FALSE}
  if (is.null(X_unval) & is.null(X_val)) {errorsX <- FALSE}
  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  # Add the B spline basis ------------------------------------------
  if (errorsX) {
    sn <- ncol(data[, Bspline])
    if(0 %in% colSums(data[c(1:n), Bspline])) {
      warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

      return(list(Coefficients = data.frame(Coefficient = NA,
                                            SE = NA),
                  converged = FALSE,
                  se_converged = NA,
                  converged_msg = "B-spline error",
                  iterations = 0))
    }
  }
  # ------------------------------------------ Add the B spline basis

  if (is.null(theta_pred)) {
    theta_pred <- c(X_val, C)
    message("Analysis model assumed main effects only.")
  }

  if (is.null(gamma_pred) & errorsY) {
    gamma_pred <- c(X_unval, Y_val, X_val, C)
    message("Outcome error model assumed main effects only.")
  }

  pred <- unique(c(theta_pred, gamma_pred))

  if (errorsX) {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, c(X_val)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- c(X_val)

    if (errorsY) {
      # Save static (X*,Y*,X,Y,C) since they don't change ---------------
      comp_dat_val <- data[c(1:n), c(Y_unval, pred, Bspline)] # c(Y_unval, X_unval, C, Bspline, X_val, Y_val)
      comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
      comp_dat_val <- comp_dat_val[, c(Y_unval, pred, Bspline, "k")]
      comp_dat_val <- data.matrix(comp_dat_val)

      # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
      # one row per x) --------------------------------------------------
      suppressWarnings(comp_dat_unval <- cbind(data[-c(1:n), c(Y_unval, setdiff(x = pred, y = c(Y_val, X_val)), Bspline)],
                                               x_obs_stacked))
      comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
      comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
      colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
      comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1),
                                          k = rep(rep(seq(1, m), each = (N - n)), times = 2)))
      comp_dat_unval <- comp_dat_unval[, c(Y_unval, pred, Bspline, "k")]

      comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
    } else {
      # Save static (X*,X,Y,C) since they don't change ---------------
      comp_dat_val <- data[c(1:n), c(pred, Bspline)]
      comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
      comp_dat_val <- comp_dat_val[, c(pred, Bspline, "k")]
      comp_dat_val <- data.matrix(comp_dat_val)

      # (m x n)xd vectors of each (one column per person, one row per x) --
      suppressWarnings(
        comp_dat_unval <- data.matrix(
          cbind(data[-c(1:n), c(setdiff(x = pred, y = c(X_val)), Bspline)],
                x_obs_stacked,
                k = rep(seq(1, m), each = (N - n)))
          )
        )
      comp_dat_unval <- comp_dat_unval[, c(pred, Bspline, "k")]

      comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
    }

    # Initialize B-spline coefficients {p_kj}  ------------
    ## Numerators sum B(Xi*) over k = 1,...,m -------------
    ## Save as p_val_num for updates ----------------------
    ## (contributions don't change) -----------------------
    p_val_num <- rowsum(x = comp_dat_val[, Bspline], group = comp_dat_val[, "k"], reorder = TRUE)
    prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))
  } else if (errorsY) {
    # Save static (Y*,X,Y,C) since they don't change ------------------
    comp_dat_val <- data.matrix(data[c(1:n), c(Y_unval, pred)])

    # Create duplicate rows of each person (one each for y = 0/1) -----
    comp_dat_unval <- data[-c(1:n), c(Y_unval, setdiff(x = pred, y = c(Y_val)))]
    comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
    comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
    colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <-
      colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
    comp_dat_unval <- data.matrix(rbind(comp_dat_y0, comp_dat_y1))

    # Stack complete data: --------------------------------------------
    ## n rows for the n subjects in Phase II (1 each) -----------------
    ## 2 * (N - n) for the (N - n) subjects from Phase I (2 each) -----
    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
  }

  theta_formula <- as.formula(paste0(Y_val, "~", paste(theta_pred, collapse = "+")))
  theta_design_mat <- cbind(int = 1, comp_dat_all[, theta_pred])

  if (errorsY) {
    gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
    gamma_design_mat <- cbind(int = 1, comp_dat_all[, gamma_pred])
  } else {
    gamma_design_mat <- matrix(1, nrow = nrow(comp_dat_all), ncol = 1)
  }

  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params <- "Zero"
  }

  if(initial_lr_params == "Zero") {
    prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
    prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)
  } else if(initial_lr_params == "Complete-data") {
    prev_theta <- theta0 <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(data[c(1:n), ]))$coefficients, ncol = 1)
    if (errorsY) {
      prev_gamma <- gamma0 <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(data[c(1:n), ]))$coefficient, ncol = 1)
    }
  }

  y0 <- comp_dat_unval[, Y_val] == 0

  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1

  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric(theta_design_mat[-c(1:n), ] %*% prev_theta)
    pY_X <- 1/(1 + exp(-mu_theta))
    pY_X[y0] <- 1-pY_X[y0]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    if (errorsY) {
      I_ystar0 <- comp_dat_unval[, Y_unval] == 0
      pYstar <- 1 / (1 + exp(- as.numeric(gamma_design_mat[-c(1:n), ] %*% prev_gamma)))
      pYstar[I_ystar0] <- 1 - pYstar[I_ystar0]

    } else {
      pYstar <- matrix(1, nrow(gamma_design_mat))
    }
    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------
    if (errorsX) {
      if (errorsY) {
        ### p_kj ------------------------------------------------------
        ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
        ### multiply by the B-spline terms
        pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ] * comp_dat_unval[, Bspline]
        ### ---------------------------------------------------------- p_kj
      } else {
        ### p_kj ----------------------------------------------------------
        ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
        ### multiply by the B-spline terms
        pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 1), ] * comp_dat_unval[, Bspline]
        ### ---------------------------------------------------------- p_kj
      }
      ### ----------------------------------------------------- P(X|X*)
    } else { pX <- rep(1, times = nrow(comp_dat_unval)) }
    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Update numerator ----------------------------------------------
    ### P(Y|X,C)*P(Y*|X*,Y,X,C)p_kjB(X*) ------------------------------
    psi_num <- pY_X * pYstar * pX
    ### Update denominator --------------------------------------------
    #### Sum up all rows per id (e.g. sum over xk/y) ------------------
    if (errorsY & errorsX) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2 * m))
    } else if (errorsX) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
    } else if (errorsY) {
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2))
    }
    #### Then sum over the sn splines ---------------------------------
    #### Same ordering as psi_t_num_sumover_k_y, just only 1 column ---
    psi_denom <- rowSums(psi_denom)
    #### Avoid NaN resulting from dividing by 0 -----------------------
    psi_denom[psi_denom == 0] <- 1
    ### And divide them! ----------------------------------------------
    psi_t <- psi_num / psi_denom
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
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
      }
    } else if (errorsY) {
      w_t <- psi_t
      ## Update the u_kji for unvalidated subjects --------------------
      ## by summing over Y = 0/1 w/i each i, k ------------------------
      ## add top half of psi_t (y = 0) to bottom half (y = 1) ---------
      u_t <- psi_t[c(1:(N - n))] + psi_t[-c(1:(N - n))]
      ## -------------------- Update the u_kji for unvalidated subjects
    }
    # ---------------------------------------------------------- E Step
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted logistic regression ----------------
    ### Gradient ------------------------------------------------------
    mu <- theta_design_mat %*% prev_theta
    w_t <- c(rep(1, n), w_t)
    gradient_theta <- matrix(data = c(colSums(w_t * c((comp_dat_all[, Y_val] - 1 + exp(-mu) / (1 + exp(- mu)))) * theta_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c((exp(-mu) / (1 + exp(-mu))) * (exp(-mu)/(1 + exp(-mu)) - 1)) * w_t * theta_design_mat
    hessian_theta <- apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)

    new_theta <- tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta,
                          error = function(err) {
                            matrix(NA, nrow = nrow(prev_theta))
                          }
    )
    if (any(is.na(new_theta))) {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta) < TOL
    ## --------------------------------------------------- Update theta
    ###################################################################
    if (errorsY) {
      ## Update gamma using weighted logistic regression ----------------
      mu <- gamma_design_mat %*% prev_gamma
      gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(-mu) / (1 + exp(-mu)))) * gamma_design_mat)), ncol = 1)
      ### ------------------------------------------------------ Gradient
      ### Hessian -------------------------------------------------------
      post_multiply <- c(w_t * (exp(- mu) / (1 + exp(- mu))) * (exp(-mu) / (1 + exp(-mu)) - 1)) * gamma_design_mat
      hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
      new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
                            error = function(err) {
                              matrix(NA, nrow = nrow(prev_gamma))
                            }
      )
      if (any(is.na(new_gamma))) {
        suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      }
      # Check for convergence -----------------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
      ## ---------------- Update gamma using weighted logistic regression
    } else {
      gamma_conv <- TRUE
    }
    ###################################################################
    if (errorsX) {
      ## Update {p_kj} --------------------------------------------------
      ### Update numerators by summing u_t over i = 1, ..., N -----------
      new_p_num <- p_val_num +
        rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
      new_p <- t(t(new_p_num) / colSums(new_p_num))
      ### Check for convergence -----------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
      ## -------------------------------------------------- Update {p_kj}
    } else {
      p_conv <- TRUE
    }
    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################

    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE

    # Update values for next iteration  -------------------------------
    it <- it + 1
    prev_theta <- new_theta
    if (errorsY) { prev_gamma <- new_gamma }
    if (errorsX) { prev_p <- new_p }
    #  ------------------------------- Update values for next iteration
  }

  if(!CONVERGED) {
    if(it > MAX_ITER) CONVERGED_MSG = "MAX_ITER reached"

    return(list(Coefficients = data.frame(Coefficient = matrix(NA, nrow = nrow(prev_theta)),
                                          SE = NA),
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = NA))
  }
  rownames(new_theta) <- c("Intercept", theta_pred)
  if (errorsY) { rownames(new_gamma) <- c("Intercept", gamma_pred) }

  if(CONVERGED) CONVERGED_MSG <- "Converged"
  # ---------------------------------------------- Estimate theta using EM
  if(noSE) {
    if (!errorsX) {
      new_p <- p_val_num <- matrix(NA, nrow = 1, ncol = 1)
    }

    if (!errorsY) {
      new_gamma <- NA
    }
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik_rw(N = N,
                                               n = n,
                                               Y_unval = Y_unval,
                                               Y_val = Y_val,
                                               X_unval = X_unval,
                                               X_val = X_val,
                                               C = C,
                                               Bspline = Bspline,
                                               comp_dat_all = comp_dat_all,
                                               theta_pred = theta_pred,
                                               gamma_pred = gamma_pred,
                                               theta = new_theta,
                                               gamma = new_gamma,
                                               p = new_p)

    return(list(Coefficients = data.frame(Coefficient = new_theta,
                                          SE = NA),
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_theta))
  } else {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_N <- h_N_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    if (!errorsX) {
      new_p <- p_val_num <- matrix(NA, nrow = 1, ncol = 1)
    }

    if (!errorsY) {
      new_gamma <- NA
    }

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik_rw(N = N,
                                            n = n,
                                            Y_unval = Y_unval,
                                            Y_val = Y_val,
                                            X_unval = X_unval,
                                            X_val = X_val,
                                            C = C,
                                            Bspline = Bspline,
                                            comp_dat_all = comp_dat_all,
                                            theta_pred = theta_pred,
                                            gamma_pred = gamma_pred,
                                            theta = new_theta,
                                            gamma = new_gamma,
                                            p = new_p)

    I_theta <- matrix(od_loglik_theta, nrow = nrow(new_theta), ncol = nrow(new_theta))

    single_pert_theta <- sapply(X = seq(1, ncol(I_theta)),
                                FUN = pl_theta_rw,
                                theta = new_theta,
                                h_N = h_N,
                                n = n,
                                N = N,
                                Y_unval = Y_unval,
                                Y_val = Y_val,
                                X_unval = X_unval,
                                X_val = X_val,
                                C = C,
                                Bspline = Bspline,
                                comp_dat_all = comp_dat_all,
                                theta_pred = theta_pred,
                                gamma_pred = gamma_pred,
                                gamma0 = new_gamma,
                                p0 = new_p,
                                p_val_num = p_val_num,
                                TOL = TOL,
                                MAX_ITER = MAX_ITER)

    if (any(is.na(single_pert_theta))) {
      I_theta <- matrix(NA, nrow = nrow(new_theta), ncol = nrow(new_theta))
      SE_CONVERGED <- FALSE
    } else {
      spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)),
                         ncol = ncol(I_theta),
                         byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta <- I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED <- TRUE
    }

    for (c in 1:ncol(I_theta)) {
      pert_theta <- new_theta
      pert_theta[c] <- pert_theta[c] + h_N
      double_pert_theta <- sapply(X = seq(c, ncol(I_theta)),
                                  FUN = pl_theta_rw,
                                  theta = pert_theta,
                                  h_N = h_N,
                                  n = n,
                                  N = N,
                                  Y_unval = Y_unval,
                                  Y_val = Y_val,
                                  X_unval = X_unval,
                                  X_val = X_val,
                                  C = C,
                                  Bspline = Bspline,
                                  comp_dat_all = comp_dat_all,
                                  theta_pred = theta_pred,
                                  gamma_pred = gamma_pred,
                                  gamma0 = new_gamma,
                                  p0 = new_p,
                                  p_val_num = p_val_num,
                                  MAX_ITER = MAX_ITER,
                                  TOL = TOL)
      dpt <- matrix(0, nrow = nrow(I_theta), ncol = ncol(I_theta))
      dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta)) {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1]
      }
      I_theta <- I_theta + dpt
    }

    I_theta <- h_N ^ (- 2) * I_theta

    cov_theta <- tryCatch(expr = - solve(I_theta),
                          error = function(err) {
                            matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta))
                          }
    )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # if(any(diag(cov_theta) < 0)) {
    #   warning("Negative variance estimate. Increase the h_N_scale parameter and repeat variance estimation.")
    #   SE_CONVERGED <- FALSE
    # }
    se_theta <- tryCatch(expr = sqrt(diag(cov_theta)),
                            warning = function(w) {
                              matrix(NA, nrow = nrow(prev_theta))
                            }
    )
    if (any(is.na(se_theta))) { SE_CONVERGED <- FALSE} else { TRUE }
    return(list(Coefficients = data.frame(Coefficient = new_theta,
                                          SE = se_theta),
                I_theta = I_theta,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_theta))
  }
}
