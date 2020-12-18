#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems
#'
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2020)
#
#'
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index)
#' @param Y_val Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param C (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Validated Columns with the validation indicator (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param data A dataframe with one row per subject containing columns: Y_unval, Y_val, X_unval, X_val, C, Validated, and Bspline.
#' @param initial_lr_params Initial values for parametric model parameters. Choices include (1) "Zero" (non-informative starting values) or (2) "Complete-data" (estimated based on validated subjects only)
#' @param h_n_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, defaults to 1.
#' @param noSE Indicator for whether standard errors are desired. Defaults to noSE = FALSE.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{Coefficients}{dataframe with final coefficient and standard error estimates (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates.}
#' \item{all_iter_theta}{dataframe containing parameter estimates from all iterations of the EM algorithm.}
#' \item{all_iter_od_loglik}{dataframe containing the value of the observed-data log-likelihood from all iterations of the EM algorithm.}
#' @export

logreg2ph <- function(Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL,
                      Validated = NULL, Bspline = NULL, data, initial_lr_params = "Zero",
                      h_n_scale = NULL, noSE = FALSE, TOL = 1E-4, MAX_ITER = 1000)
{
  rescale <- FALSE

  N <- nrow(data)
  n <- sum(data[, Validated])

  if(is.null(h_n_scale)) {
    h_n_scale <- 1
  }

  # Reorder so that the n validated subjects are first ------------
  data <- data[order(as.numeric(data[, Validated]), decreasing = TRUE), ]

  # Add the B spline basis ------------------------------------------
  sn <- ncol(data[, Bspline])
  if(0 %in% colSums(data[c(1:n), Bspline])) {
    warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

    return(list(Coefficients = data.frame(Coefficient = NA,
                                          SE = NA),
                converged = FALSE,
                se_converged = FALSE,
                converged_msg = "B-spline error",
                iterations = 0,
                all_iter_theta = NA,
                all_iter_od_loglik = NA,
                od_loglik_at_conv = NA))
  }
  # ------------------------------------------ Add the B spline basis

  # Standardize X_val, X_unval and C to N(0,1) ----------------------
  re_data <- data
  if(rescale) {
    ## Shift by the sample mean ---------------------------------------
    re_shift <- c(0, as.numeric(colMeans(data[, c(X_val, X_unval, C)], na.rm = TRUE)))
    ## Scale inversely by sample standard deviation -------------------
    re_scale <- c(1, as.numeric(apply(data[, c(X_val, X_unval, C)], MARGIN = 2, FUN = sd, na.rm = TRUE)))
    for (p in 1:length(c(X_val, X_unval, C)))
    {
      re_data[,c(X_val, X_unval, C)[p]] <- (re_data[,c(X_val, X_unval, C)[p]] - re_shift[p + 1]) / re_scale[p + 1]
    }
  }
  # ---------------------- Standardize X_val, X_unval and C to N(0,1)

  # Save distinct X -------------------------------------------------
  x_obs <- data.frame(unique(re_data[1:n, c(X_val)]))
  x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
  m <- nrow(x_obs)
  x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
  x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
  colnames(x_obs) <- colnames(x_obs_stacked) <- c(X_val)

  # Save static (X*,Y*,X,Y,C) since they don't change ---------------
  comp_dat_val <- re_data[c(1:n), c(Y_unval, X_unval, C, Bspline, X_val, Y_val)]
  comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
  comp_dat_val <- comp_dat_val[, c(Y_unval, X_unval, C, Bspline, X_val, Y_val, "k")]
  comp_dat_val <- data.matrix(comp_dat_val)
  # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
  # one row per x) --------------------------------------------------
  suppressWarnings(comp_dat_unval <- cbind(re_data[-c(1:n), c(Y_unval, X_unval, C, Bspline)],
                                           x_obs_stacked))
  comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
  comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
  colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
  comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1),
                                      k = rep(rep(seq(1, m), each = (N - n)), times = 2)))

  comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)

  theta_formula <- as.formula(paste0(Y_val, "~", paste(c(X_val, C), collapse = "+")))
  theta_design_mat <- cbind(int = 1, comp_dat_all[, c(X_val, C)])

  gamma_formula <- as.formula(paste0(Y_unval, "~", paste(c(X_unval, Y_val, X_val, C), collapse = "+")))
  gamma_design_mat <- cbind(int = 1, comp_dat_all[, c(X_unval, Y_val, X_val, C)])

  # Standardize Y_val in gamma_design_mat to N(0,1) -----------------
  ## Shift by the sample mean ---------------------------------------
  ## Scale inversely by sample standard deviation -------------------
  ## Create artificially scaled data set ----------------------------
  if(rescale) {
    gamma_design_mat[, Y_val] <- (gamma_design_mat[, Y_val] - mean(gamma_design_mat[, Y_val]))/sd(gamma_design_mat[, Y_val])
  }
  ## Everything else was already standardized above -----------------
  # -------------------------------  Standardize Y_val in gamma_design_mat to N(0,1)

  # If unvalidated variable was left blank, assume error-free -------
  ## Need to write simplification here

  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data"))) {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params <- "Zero"
  }
  if(initial_lr_params == "Zero") {
    num_pred <- length(X_val) + length(C) #preds in analysis model --
    prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
    prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)
  }
  if(initial_lr_params == "Complete-data") {
    prev_theta <- theta0 <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(re_data[c(1:n), ]))$coefficients, ncol = 1)
    prev_gamma <- gamma0 <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(re_data[c(1:n), ]))$coefficient, ncol = 1)
  }

  # Initialize B-spline coefficients {p_kj}  ------------
  ## Numerators sum B(Xi*) over k = 1,...,m -------------
  ## Save as p_val_num for updates ----------------------
  ## (contributions don't change) -----------------------
  p_val_num <- rowsum(x = comp_dat_val[, Bspline], group = comp_dat_val[, "k"], reorder = TRUE)
  prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))

  y0 <- comp_dat_unval[, Y_val] == 0
  ystar0 <- comp_dat_unval[, Y_unval] == 0

  CONVERGED <- FALSE
  CONVERGED_MSG <- "Unknown"
  it <- 1

  all_theta <- all_od_loglik <- vector()

  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED) {
    # Save values
    all_theta <- append(all_theta, prev_theta)
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric(theta_design_mat[-c(1:n), ] %*% prev_theta)
    pY_X <- 1/(1 + exp(-mu_theta))
    pY_X[y0] <- 1-pY_X[y0]
    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    pYstar <- 1/(1 + exp(-as.numeric(gamma_design_mat[-c(1:n), ] %*% prev_gamma)))
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
    #### Sum up all rows per id (e.g. sum over xk/y) ------------------
    #### reorder = TRUE returns them in ascending order of i ----------
    #### (rather than in order of encounter) --------------------------
    psi_t_num_sumover_k_y <- rowsum(psi_t_num, group = rep(seq(1, (N - n)), times = 2 * m))
    #### Then sum over the sn splines ---------------------------------
    #### Same ordering as psi_t_num_sumover_k_y, just only 1 column ---
    psi_t_num_sumover_k_y_j <- rowSums(psi_t_num_sumover_k_y)
    #### Avoid NaN resulting from dividing by 0 -----------------------
    psi_t_num_sumover_k_y_j[psi_t_num_sumover_k_y_j == 0] <- 1
    ### And divide them! ----------------------------------------------
    psi_t <- psi_t_num / psi_t_num_sumover_k_y_j
    ## ------------------- Update the psi_kyji for unvalidated subjects
    ###################################################################
    ## Update the w_kyi for unvalidated subjects ----------------------
    ## by summing across the splines/ columns of psi_t ----------------
    ## w_t is ordered by i = 1, ..., N --------------------------------
    w_t <- rowSums(psi_t)
    ## ---------------------- Update the w_kyi for unvalidated subjects
    ## For validated subjects, w_t = I(Xi=xk) so make them all 0 ------
    #w_t[rep(data[,Validated], 2*m)] <- 0
    ## then place a 1 in the w_t_val positions ------------------------
    #w_t[w_t_val] <- 1
    ## Check: w_t sums to 1 over within i -----------------------------
    # table(rowSums(rowsum(w_t, group = rep(rep(seq(1,N), times = m), times = 2))))

    ## Update the u_kji for unvalidated subjects ----------------------
    ## by summing over Y = 0/1 w/i each i, k --------------------------
    ## add top half of psi_t (y = 0) to bottom half (y = 1) -----------
    u_t <- psi_t[c(1:(m * (N - n))),] + psi_t[-c(1:(m * (N - n))), ]
    ## make u_t for the (1:n) validated subjects = 0 ----------------
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
    w_t <- c(rep(1,n), w_t)
    gradient_theta <- matrix(data = c(colSums(w_t * c((comp_dat_all[, Y_val]-1 + exp(-mu) / (1+exp(-mu)))) * theta_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c((exp(-mu) / (1 + exp(-mu))) * (exp(-mu)/(1 + exp(-mu)) - 1)) * w_t * theta_design_mat
    hessian_theta <- apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)

    new_theta <- tryCatch(expr = prev_theta - solve(hessian_theta) %*% gradient_theta,
                          error = function(err) {
                            matrix(NA, nrow = nrow(prev_theta))
                          }
    )
    if (TRUE %in% is.na(new_theta)) {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta) < TOL
    ## --------------------------------------------------- Update theta
    ###################################################################
    ## Update gamma using weighted logistic regression ----------------
    mu <- gamma_design_mat %*% prev_gamma
    gradient_gamma <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(-mu) / (1 + exp(-mu)))) * gamma_design_mat)), ncol = 1)
    ### ------------------------------------------------------ Gradient
    ### Hessian -------------------------------------------------------
    post_multiply <- c(w_t*(exp(-mu)/(1+exp(-mu))) * (exp(-mu) / (1 + exp(-mu)) - 1)) * gamma_design_mat
    hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    new_gamma <- tryCatch(expr = prev_gamma - solve(hessian_gamma) %*% gradient_gamma,
                          error = function(err) {
                            matrix(NA, nrow = nrow(prev_gamma))
                          }
    )
    if (TRUE %in% is.na(new_gamma)) {
      suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
    }
    # Check for convergence -----------------------------------------
    gamma_conv <- abs(new_gamma - prev_gamma) < TOL
    ## ---------------- Update gamma using weighted logistic regression
    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    ### Update numerators by summing u_t over i = 1, ..., N -----------
    new_p_num <- p_val_num +
      rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
    new_p <- t(t(new_p_num) / colSums(new_p_num))
    ### Check for convergence -----------------------------------------
    p_conv <- abs(new_p - prev_p) < TOL
    ## -------------------------------------------------- Update {p_kj}
    ###################################################################
    # M Step ----------------------------------------------------------

    all_od_loglik <- append(all_od_loglik,
                            observed_data_loglik(N = N, n = n,
                                                 Y_unval=Y_unval, Y_val=Y_val,
                                                 X_unval=X_unval, X_val=X_val,
                                                 C=C, Bspline=Bspline,
                                                 comp_dat_all = comp_dat_all,
                                                 theta = new_theta,
                                                 gamma = new_gamma,
                                                 p = new_p))

    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1) CONVERGED <- TRUE

    it <- it + 1

    # Update values for next iteration  -------------------------------
    prev_theta <- new_theta
    prev_gamma <- new_gamma
    prev_p <- new_p
  }

  if(!CONVERGED) {
    if(it > MAX_ITER) CONVERGED_MSG = "MAX_ITER reached"

    return(list(Coefficients = data.frame(Coefficient = matrix(NA, nrow = nrow(prev_theta)),
                                          SE = NA),
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                all_iter_theta = all_theta,
                all_iter_od_loglik = all_od_loglik,
                od_loglik_at_conv = NA))
  }
  if(CONVERGED) CONVERGED_MSG <- "Converged"
  # ---------------------------------------------- Estimate theta using EM
  if(noSE) {
    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N, n = n,
                                            Y_unval=Y_unval, Y_val=Y_val,
                                            X_unval=X_unval, X_val=X_val,
                                            C=C, Bspline=Bspline,
                                            comp_dat_all = comp_dat_all,
                                            theta = new_theta,
                                            gamma = new_gamma,
                                            p = new_p)

    rownames(new_theta) <- c("Intercept", colnames(theta_design_mat)[-1])
    re_theta <- new_theta
    rownames(new_gamma) <- c("Intercept", X_unval, Y_val, X_val, C)
    re_gamma <- new_gamma
    if(rescale) {
      re_theta[c(2:(1 + length(c(X_val, C))))] <- re_theta[c(2:(1 + length(c(X_val, C))))] / re_scale[c(2:(1 + length(c(X_val, C))))]
      re_theta[1] <- re_theta[1] - sum(re_theta[c(2:(1 + length(c(X_val, C))))] * re_shift[c(2:(1 + length(c(X_val, C))))])

      re_gamma[c(2:(1 + length(c(X_val, C))))] <- re_gamma[c(2:(1 + length(c(X_val, C))))] / re_scale[c(2:(1 + length(c(X_val, C))))]
      re_gamma[1] <- re_gamma[1] - sum(re_gamma[c(2:(1 + length(c(X_val, C))))] * re_shift[c(2:(1 + length(c(X_val, C))))])
    }
    return(list(Coefficients = data.frame(Coefficient = re_theta,
                                          SE = NA),
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                all_iter_theta = all_theta,
                all_iter_od_loglik = all_od_loglik,
                od_loglik_at_conv = od_loglik_theta))
  } else {
    rownames(new_theta) <- c("Intercept", colnames(theta_design_mat)[-1])
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_n <- h_n_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N, n = n,
                                            Y_unval=Y_unval, Y_val=Y_val,
                                            X_unval=X_unval, X_val=X_val,
                                            C=C, Bspline=Bspline,
                                            comp_dat_all = comp_dat_all,
                                            theta = new_theta,
                                            gamma = new_gamma,
                                            p = new_p)

    I_theta <- matrix(od_loglik_theta, nrow = nrow(new_theta), ncol = nrow(new_theta))

    #single perts: 271.658 sec elapsed
    single_pert_theta <- sapply(X = seq(1,ncol(I_theta)),
                                FUN = pl_theta,
                                theta = new_theta, h_n = h_n, n = n, N = N,
                                Y_unval=Y_unval, Y_val=Y_val, X_unval=X_unval, X_val=X_val,
                                C=C, Bspline=Bspline,
                                comp_dat_all = comp_dat_all, gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num,
                                TOL = TOL, MAX_ITER = MAX_ITER)

    if (NA %in% single_pert_theta) {
      I_theta <- matrix(NA, nrow = nrow(new_theta), ncol = nrow(new_theta))
      SE_CONVERGED <- FALSE
    } else {
      spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)), ncol = ncol(I_theta), byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta <- I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED <- TRUE
    }

    for (c in 1:ncol(I_theta)) {
      pert_theta <- new_theta; pert_theta[c] <- pert_theta[c] + h_n
      double_pert_theta <- sapply(X = seq(c,ncol(I_theta)), FUN = pl_theta, theta = pert_theta, h_n = h_n, n = n, N = N,
                                  Y_unval=Y_unval, Y_val=Y_val, X_unval=X_unval, X_val=X_val, C=C, Bspline=Bspline,
                                  comp_dat_all = comp_dat_all, gamma0 = new_gamma, p0 = new_p, p_val_num = p_val_num, MAX_ITER = MAX_ITER, TOL = TOL)
      dpt <- matrix(0, nrow = nrow(I_theta), ncol = ncol(I_theta))
      dpt[c,c] <- double_pert_theta[1] #Put double on the diagonal
      if(c < ncol(I_theta)) {
        dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1] #And fill the others in on the cth row/ column
      }
      I_theta <- I_theta + dpt
    }

    I_theta <- h_n ^ ( - 2) * I_theta

    cov_theta <- tryCatch(expr = - solve(I_theta),
                          error = function(err) {
                            matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta))
                          }
    )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # Scale everything back ------------------------------------------------
    re_theta <- new_theta
    if(TRUE %in% (diag(cov_theta) < 0)) {
      warning("Negative variance estimate. Increase the h_n_scale parameter and repeat variance estimation.")
      SE_CONVERGED <- FALSE
    }
    re_se_theta <- tryCatch(expr = sqrt(diag(cov_theta)),
                            warning = function(w) {
                              matrix(NA, nrow = nrow(prev_theta))
                            }
    )
    re_gamma <- new_gamma

    if(rescale) {
      re_theta[c(2:(1+length(c(X_val, C))))] <- re_theta[c(2:(1 + length(c(X_val, C))))] / re_scale[c(2:(1 + length(c(X_val, C))))]
      re_theta[1] <- re_theta[1] - sum(re_theta[c(2:(1 + length(c(X_val, C))))] * re_shift[c(2:(1 + length(c(X_val, C))))])

      re_se_theta[c(2:(1 + length(c(X_val, C))))] <- re_se_theta[c(2:(1 + length(c(X_val, C))))] / re_scale[c(2:(1 + length(c(X_val, C))))]
      re_se_theta[1] <- cov_theta[1,1] + sum(diag(cov_theta)[c(2:(1 + length(c(X_val, C))))] * (re_shift[c(2:(1 + length(c(X_val, C))))] / re_scale[c(2:(1 + length(c(X_val, C))))]) ^ 2)
      for (p1 in 1:ncol(cov_theta)) {
        for (p2 in p1:ncol(cov_theta)) {
          if(p1 < p2 & p1 == 1) {
            re_se_theta[1] <- re_se_theta[1] - (re_shift[p2] / re_scale[p2]) * cov_theta[p1, p2]
          }
          if(p1 < p2 & p1 > 1) {
            re_se_theta[1] <- re_se_theta[1] + (re_shift[p1] / re_scale[p1]) * (re_shift[p2] / re_scale[p2]) * cov_theta[p1, p2]
          }
        }
      }
      re_se_theta[1] <- sqrt(re_se_theta[1])

      re_gamma[c(2:(1 + length(c(X_val, C))))] <- re_gamma[c(2:(1 + length(c(X_val, C))))] / re_scale[c(2:(1 + length(c(X_val, C))))]
      re_gamma[1] <- re_gamma[1] - sum(re_gamma[c(2:(1 + length(c(X_val, C))))] * re_shift[c(2:(1 + length(c(X_val, C))))])
    }
    # ------------------------------------------------ Scale everything back
    return(list(Coefficients = data.frame(Coefficient = re_theta,
                                          SE = re_se_theta),
                I_theta = I_theta,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                all_iter_theta = all_theta,
                all_iter_od_loglik = all_od_loglik,
                od_loglik_at_conv = od_loglik_theta))
  }
}
