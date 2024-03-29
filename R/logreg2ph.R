#' Sieve maximum likelihood estimator (SMLE) for two-phase logistic regression problems
#' This function returns the sieve maximum likelihood estimators (SMLE) for the logistic regression model from Lotspeich et al. (2021)
#'
#' @param Y_unval Column name with the unvalidated outcome. If \code{Y_unval} is null, the outcome is assumed to be error-free.
#' @param Y_val Column name with the validated outcome.
#' @param X_unval Column name(s) with the unvalidated predictors.  If \code{X_unval} and \code{X_val} are \code{null}, all precictors are assumed to be error-free.
#' @param X_val Column name(s) with the validated predictors. If \code{X_unval} and \code{X_val} are \code{null}, all precictors are assumed to be error-free.
#' @param C (Optional) Column name(s) with additional error-free covariates.
#' @param Validated Column name with the validation indicator. The validation indicator can be defined as \code{Validated = 1} or \code{TRUE} if the subject was validated and \code{Validated = 0} or \code{FALSE} otherwise.
#' @param Bspline Vector of column names containing the B-spline basis functions.
#' @param data A dataframe with one row per subject containing columns: \code{Y_unval}, \code{Y_val}, \code{X_unval}, \code{X_val}, \code{C}, \code{Validated}, and \code{Bspline}.
#' @param theta_pred Vector of columns in \code{data} that pertain to the predictors in the analysis model.
#' @param gamma_pred Vector of columns in \code{data} that pertain to the predictors in the outcome error model.
#' @param initial_lr_params Initial values for parametric model parameters. Choices include (1) \code{"Zero"} (non-informative starting values) or (2) \code{"Complete-data"} (estimated based on validated subjects only)
#' @param h_N_scale Size of the perturbation used in estimating the standard errors via profile likelihood. If none is supplied, default is `h_N_scale = 1`.
#' @param noSE Indicator for whether standard errors are desired. Defaults to \code{noSE = FALSE}.
#' @param TOL Tolerance between iterations in the EM algorithm used to define convergence.
#' @param MAX_ITER Maximum number of iterations allowed in the EM algorithm.
#' @return
#' \item{coeff}{dataframe with final coefficient and standard error estimates (where applicable) for the analysis model.}
#' \item{outcome_err_coeff}{dataframe with final coefficient estimates for the outcome error model.}
#' \item{Bspline_coeff}{dataframe with final B-spline coefficient estimates (where applicable).}
#' \item{vcov}{variance-covarianced matrix for \code{coeff} (where applicable).}
#' \item{converged}{indicator of EM algorithm convergence for parameter estimates.}
#' \item{se_converged}{indicator of standard error estimate convergence.}
#' \item{converged_msg}{(where applicable) description of non-convergence.}
#' \item{iterations}{number of iterations completed by EM algorithm to find parameter estimates.}
#' \item{od_loglik_at_conv}{value of the observed-data log-likelihood at convergence.}
#' @export

logreg2ph <- function(Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL,
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
  if (is.null(Y_unval))
  {errorsY <- FALSE}

  if (is.null(X_unval) & is.null(X_val))
  {errorsX <- FALSE}

  ## ------ If unvalidated variable was left blank, assume error-free
  # ----------------------------------------- Determine error setting

  # Add the B spline basis ------------------------------------------
  if (errorsX)
  {
    sn <- ncol(data[, Bspline])
    if(0 %in% colSums(data[c(1:n), Bspline]))
    {
      warning("Empty sieve in validated data. Reconstruct B-spline basis and try again.", call. = FALSE)

      return(list(coeff = data.frame(coeff = NA, se = NA),
                  outcome_err_coeff = data.frame(coeff = NA, se = NA),
                  Bspline_coeff = NA,
                  vcov = NA,
                  converged = NA,
                  se_converged = NA,
                  converged_msg = "B-spline error",
                  iterations = 0,
                  od_loglik_at_conv = NA))
    }

  }

  # ------------------------------------------ Add the B spline basis

  if (is.null(theta_pred))
  {
    theta_pred <- c(X_val, C)
    message("Analysis model assumed main effects only.")
  }


  if (is.null(gamma_pred) & errorsY)
  {
    gamma_pred <- c(X_unval, Y_val, X_val, C)
    message("Outcome error model assumed main effects only.")
  }


  pred <- unique(c(theta_pred, gamma_pred))

  if (errorsX & errorsY)
  {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, c(X_val)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- c(X_val)

    # Save static (X*,Y*,X,Y,C) since they don't change ---------------
    comp_dat_val <- data[c(1:n), c(Y_unval, X_unval, C, Bspline, X_val, Y_val)]
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

    # Initialize B-spline coefficients {p_kj}  ------------
    ## Numerators sum B(Xi*) over k = 1,...,m -------------
    ## Save as p_val_num for updates ----------------------
    ## (contributions don't change) -----------------------
    p_val_num <- rowsum(x = comp_dat_val[, Bspline], group = comp_dat_val[, "k"], reorder = TRUE)
    prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))
  }
  else if (errorsX)
  {
    # Save distinct X -------------------------------------------------
    x_obs <- data.frame(unique(data[1:n, c(X_val)]))
    x_obs <- data.frame(x_obs[order(x_obs[, 1]), ])
    m <- nrow(x_obs)
    x_obs_stacked <- do.call(rbind, replicate(n = (N - n), expr = x_obs, simplify = FALSE))
    x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[, 1]), ])
    colnames(x_obs) <- colnames(x_obs_stacked) <- c(X_val)

    # Save static (X*,X,Y,C) since they don't change ---------------
    comp_dat_val <- data[c(1:n), c(Y_val, pred, Bspline)]
    comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
    comp_dat_val <- comp_dat_val[, c(Y_val, pred, Bspline, "k")]
    comp_dat_val <- data.matrix(comp_dat_val)

    # (m x n)xd vectors of each (one column per person, one row per x) --
    suppressWarnings(
      comp_dat_unval <- data.matrix(
        cbind(data[-c(1:n), c(Y_val, setdiff(x = pred, y = c(X_val)), Bspline)],
          x_obs_stacked,
          k = rep(seq(1, m), each = (N - n)))
        )
      )
    comp_dat_unval <- comp_dat_unval[, c(Y_val, pred, Bspline, "k")]

    comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)

    # Initialize B-spline coefficients {p_kj}  ------------
    ## Numerators sum B(Xi*) over k = 1,...,m -------------
    ## Save as p_val_num for updates ----------------------
    ## (contributions don't change) -----------------------
    p_val_num <- rowsum(x = comp_dat_val[, Bspline], group = comp_dat_val[, "k"], reorder = TRUE)
    prev_p <- p0 <-  t(t(p_val_num) / colSums(p_val_num))
  }
  else if (errorsY)
  {
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

  if (errorsY)
  {
    gamma_formula <- as.formula(paste0(Y_unval, "~", paste(gamma_pred, collapse = "+")))
    gamma_design_mat <- cbind(int = 1, comp_dat_all[, gamma_pred])
  }


  # Initialize parameter values -------------------------------------
  ## theta, gamma ---------------------------------------------------
  if(!(initial_lr_params %in% c("Zero", "Complete-data")))
  {
    message("Invalid starting values provided. Non-informative zeros assumed.")
    initial_lr_params <- "Zero"
  }


  if(initial_lr_params == "Zero")
  {
    prev_theta <- theta0 <- matrix(0, nrow = ncol(theta_design_mat), ncol = 1)
    if (errorsY)
    {
      prev_gamma <- gamma0 <- matrix(0, nrow = ncol(gamma_design_mat), ncol = 1)
    }

    } else if(initial_lr_params == "Complete-data")
    {
      prev_theta <- theta0 <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(data[c(1:n), ]))$coefficients, ncol = 1)
      if (errorsY)
      {
        prev_gamma <- gamma0 <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(data[c(1:n), ]))$coefficient, ncol = 1)
      }

    }


    CONVERGED <- FALSE
    CONVERGED_MSG <- "Unknown"
    it <- 1

  # pre-allocate memory for loop variables
  mus_theta <- vector("numeric", nrow(theta_design_mat) * ncol(prev_theta))
  mus_gamma <- vector("numeric", nrow(gamma_design_mat) * ncol(prev_gamma))

  # Estimate theta using EM -------------------------------------------
  while(it <= MAX_ITER & !CONVERGED)
  {
    # E Step ----------------------------------------------------------
    ## Update the psi_kyji for unvalidated subjects -------------------
    ### P(Y|X) --------------------------------------------------------
    mu_theta <- as.numeric((theta_design_mat[-c(1:n), ] %*% prev_theta))
    pY_X <- 1 / (1 + exp(- mu_theta))
    I_y0 <- comp_dat_unval[, Y_val] == 0
    pY_X[I_y0] <- 1 - pY_X[I_y0]

    ### -------------------------------------------------------- P(Y|X)
    ###################################################################
    ### P(Y*|X*,Y,X) --------------------------------------------------
    if (errorsY)
    {
      pYstar <- 1 / (1 + exp(- as.numeric((gamma_design_mat[-c(1:n), ] %*% prev_gamma))))
      I_ystar0 <- comp_dat_unval[, Y_unval] == 0
      pYstar[I_ystar0] <- 1 - pYstar[I_ystar0]

    } #else {
      #pYstar <- matrix(1, nrow = nrow(comp_dat_unval))
    #}

    ### -------------------------------------------------- P(Y*|X*,Y,X)
    ###################################################################
    ### P(X|X*) -------------------------------------------------------
    # cpppX <- pXCalc(n, comp_dat_unval, errorsX, errorsY, prev_p, rep(seq(1, m), each = (N - n))-1, match(Bspline, colnames(comp_dat_unval))-1, seq(1, nrow(comp_dat_unval))-1)
    if (errorsX & errorsY)
    {
      ### p_kj ------------------------------------------------------
      ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
      ### multiply by the B-spline terms
      pX <- prev_p[rep(rep(seq(1, m), each = (N - n)), times = 2), ] * comp_dat_unval[, Bspline]
      ### ---------------------------------------------------------- p_kj
    }
    else if (errorsX)
    {
      ### p_kj ----------------------------------------------------------
      ### need to reorder pX so that it's x1, ..., x1, ...., xm, ..., xm-
      ### multiply by the B-spline terms
      pX <- prev_p[rep(seq(1, m), each = (N - n)), ]  * comp_dat_unval[, Bspline]
      ### ---------------------------------------------------------- p_kj
    } #else if (errorsY) {
      #pX <- rep(1, times = nrow(comp_dat_unval))
    #}

    # tic("cond exp cpp")
    # condExp <- conditionalExpectations(errorsX, errorsY, pX, pY_X, pYstar, N-n, m)
    # CPPw_t <- condExp[["w_t"]]
    # CPPu_t <- condExp[["u_t"]]
    # CPPpsi_t <- condExp[["psi_t"]]
    # toc()

    # R is faster than cpp

    ### ------------------------------------------------------- P(X|X*)
    ###################################################################
    ### Estimate conditional expectations -----------------------------
    if (errorsY & errorsX)
    {
      ### P(Y|X,C)P(Y*|X*,Y,X,C)p_kjB(X*) -----------------------------
      psi_num <- c(pY_X * pYstar) * pX
      ### Update denominator ------------------------------------------
      #### Sum up all rows per id (e.g. sum over xk/y) ----------------
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2 * m))
      #### Then sum over the sn splines -------------------------------
      psi_denom <- rowSums(psi_denom)
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] <- 1
      ### And divide them! --------------------------------------------
      psi_t <- psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      ### by summing across the splines/ columns of psi_t -------------
      w_t <- rowSums(psi_t)
      ### Update the u_kji for unvalidated subjects ------------------
      ### by summing over Y = 0/1 w/i each i, k ----------------------
      ### add top half of psi_t (y = 0) to bottom half (y = 1) -------
      u_t <- psi_t[c(1:(m * (N - n))), ] + psi_t[- c(1:(m * (N - n))), ]
      } else if (errorsX)
      {
      ### P(Y|X,C)p_kjB(X*) -------------------------------------------
      psi_num <- c(pY_X) * pX
      ### Update denominator ------------------------------------------
      #### Sum up all rows per id (e.g. sum over xk) ------------------
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = m))
      #### Then sum over the sn splines -------------------------------
      psi_denom <- rowSums(psi_denom)
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] <- 1
      ### And divide them! --------------------------------------------
      psi_t <- psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      ### by summing across the splines/ columns of psi_t -------------
      w_t <- rowSums(psi_t)
      } else if (errorsY)
      {
      ### P(Y|X,C)P(Y*|Y,X,C) -----------------------------------------
      #### Sum up all rows per id (e.g. sum over y) -------------------
      psi_num <- matrix(c(pY_X * pYstar), ncol = 1)
      psi_denom <- rowsum(psi_num, group = rep(seq(1, (N - n)), times = 2))
      #### Avoid NaN resulting from dividing by 0 ---------------------
      psi_denom[psi_denom == 0] <- 1
      ### And divide them! --------------------------------------------
      psi_t <- psi_num / psi_denom
      ### Update the w_kyi for unvalidated subjects -------------------
      w_t <- psi_t
    }
    # if (max(c(abs(CPPu_t - u_t), abs(CPPw_t - w_t), abs(CPPpsi_t - psi_t))) > 1e-10)
    # {
    #   warning("cpp and R are significantly different!")
    #   browser()
    # }
    ### ----------------------------- Estimate conditional expectations
    # ---------------------------------------------------------- E Step
    ###################################################################

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################
    ## Update theta using weighted logistic regression ----------------
    ### Gradient ------------------------------------------------------
    w_t <- lengthenWT(w_t, n)

    # calculateMu returns exp(-mu) / (1 + exp(-mu))
    muVector <- calculateMu(theta_design_mat, prev_theta)
    gradient_theta <- calculateGradient(w_t, n, theta_design_mat, comp_dat_all[, Y_val], muVector)
    hessian_theta <- calculateHessian(theta_design_mat, w_t, muVector, n, mus_theta);


    # tic("m-step R")
    # mu <- theta_design_mat %*% prev_theta
    # w_t <- c(rep(1, n), w_t)
    # gradient_theta_R <- matrix(data = c(colSums(w_t * c((comp_dat_all[, Y_val] - 1 + exp(-mu) / (1 + exp(- mu)))) * theta_design_mat)), ncol = 1)
    # toc()
    #
    #
    # ### ------------------------------------------------------ Gradient
    # ### Hessian -------------------------------------------------------
    # post_multiply <- c((muVector) * (muVector - 1)) * w_t * theta_design_mat
    # hessian_theta_R <- apply(theta_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
    #
    # if (max(abs(hessian_theta - hessian_theta_R)) > 1e-10)
    # {
    #   warning("cpp and R are significantly different!")
    #   browser()
    # }
    new_theta <- tryCatch(expr = prev_theta - (solve(hessian_theta) %*% gradient_theta),
      error = function(err)
      {
        matrix(NA, nrow = nrow(prev_theta))
        })
    if (any(is.na(new_theta)))
    {
      suppressWarnings(new_theta <- matrix(glm(formula = theta_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      # browser()
    }

    ### Check for convergence -----------------------------------------
    theta_conv <- abs(new_theta - prev_theta) < TOL

    ## --------------------------------------------------- Update theta
    ###################################################################
    if (errorsY)
    {
      # w_t is already the proper size

      ## Update gamma using weighted logistic regression ----------------
      muVector <- calculateMu(gamma_design_mat, prev_gamma)
      gradient_gamma <- calculateGradient(w_t, n, gamma_design_mat, comp_dat_all[, c(Y_unval)], muVector)
      hessian_gamma <- calculateHessian(gamma_design_mat, w_t, muVector, n, mus_gamma)

      # mu <- gamma_design_mat %*% prev_gamma
      # gradient_gamma_R <- matrix(data = c(colSums(w_t * c((comp_dat_all[, c(Y_unval)] - 1 + exp(- mu) / (1 + exp(- mu)))) * gamma_design_mat)), ncol = 1)

      # if (max(abs(gradient_gamma - gradient_gamma_R)) > 1e-10)
      # {
      #   warning("cpp and R are significantly different! gamma")
      #   browser()
      # }

      # ### ------------------------------------------------------ Gradient
      # ### Hessian -------------------------------------------------------
      # post_multiply <- c(w_t * muVector * (muVector - 1)) * gamma_design_mat
      # hessian_gamma <- apply(gamma_design_mat, MARGIN = 2, FUN = hessian_row, pm = post_multiply)
      new_gamma <- tryCatch(expr = prev_gamma - (solve(hessian_gamma) %*% gradient_gamma),
        error = function(err)
        {
          matrix(NA, nrow = nrow(prev_gamma))
          })
      if (any(is.na(new_gamma)))
      {
        suppressWarnings(new_gamma <- matrix(glm(formula = gamma_formula, family = "binomial", data = data.frame(comp_dat_all), weights = w_t)$coefficients, ncol = 1))
      }

      # Check for convergence -----------------------------------------
      gamma_conv <- abs(new_gamma - prev_gamma) < TOL
      ## ---------------- Update gamma using weighted logistic regression
    }
    else
    { gamma_conv <- TRUE }

    ###################################################################
    ## Update {p_kj} --------------------------------------------------
    if (errorsX & errorsY)
    {
      ### Update numerators by summing u_t over i = 1, ..., N ---------
      new_p_num <- p_val_num +
      rowsum(u_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
      new_p <- t(t(new_p_num) / colSums(new_p_num))
      ### Check for convergence ---------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
    }
    else if (errorsX)
    {
      ### Update numerators by summing u_t over i = 1, ..., N ---------
      new_p_num <- p_val_num +
      rowsum(psi_t, group = rep(seq(1, m), each = (N - n)), reorder = TRUE)
      new_p <- t(t(new_p_num) / colSums(new_p_num))
      ### Check for convergence ---------------------------------------
      p_conv <- abs(new_p - prev_p) < TOL
    }
    else
    { p_conv <- TRUE }

    ## -------------------------------------------------- Update {p_kj}

    ###################################################################
    # M Step ----------------------------------------------------------
    ###################################################################

    all_conv <- c(theta_conv, gamma_conv, p_conv)
    if (mean(all_conv) == 1)
    { CONVERGED <- TRUE }


    # Update values for next iteration  -------------------------------
    it <- it + 1
    prev_theta <- new_theta
    if (errorsY)
    { prev_gamma <- new_gamma }

    if (errorsX)
    { prev_p <- new_p }

    #  ------------------------------- Update values for next iteration
  }


  rownames(new_theta) <- c("Intercept", theta_pred)
  if (errorsY)
  {
    rownames(new_gamma) <- c("Intercept", gamma_pred)
  }


  if(!CONVERGED)
  {
    if(it > MAX_ITER)
    {
      CONVERGED_MSG = "MAX_ITER reached"
    }


    return(list(coeff = data.frame(coeff = NA, se = NA),  
                outcome_err_coeff = data.frame(coeff = NA, se = NA),  
                Bspline_coeff = NA, 
                vcov = NA,  
                converged = FALSE,  
                se_converged = NA,  
                converged_msg = "MAX_ITER reached", 
                iterations = it,  
                od_loglik_at_conv = NA))
  }

  if(CONVERGED)
  { CONVERGED_MSG <- "Converged" }

  # ---------------------------------------------- Estimate theta using EM
  if(noSE)
  {
    if (!errorsX)
    {
      new_p <- p_val_num <- matrix(NA, nrow = 1, ncol = 1)
    }

    if (!errorsY)
    {
      new_gamma <- NA
    }

    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N,
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

    return(list(coeff = data.frame(coeff = new_theta, se = NA),
                outcome_err_coeff = data.frame(coeff = new_gamma, se = NA),
                Bspline_coeff = cbind(k = comp_dat_val[, "k"], new_p),
                vcov = NA,
                converged = CONVERGED,
                se_converged = NA,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_theta))
  }
  else
  {
    # Estimate Cov(theta) using profile likelihood -------------------------
    h_N <- h_N_scale * N ^ ( - 1 / 2) # perturbation ----------------------------

    if (!errorsX)
    {
      new_p <- p_val_num <- matrix(NA, nrow = 1, ncol = 1)
    }

    if (!errorsY)
    {
      new_gamma <- NA
    }


    ## Calculate pl(theta) -------------------------------------------------
    od_loglik_theta <- observed_data_loglik(N = N,
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
      FUN = pl_theta,
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

    if (any(is.na(single_pert_theta)))
    {
      I_theta <- matrix(NA, nrow = nrow(new_theta), ncol = nrow(new_theta))
      SE_CONVERGED <- FALSE
    }
    else
    {
      spt_wide <- matrix(rep(c(single_pert_theta), times = ncol(I_theta)),
       ncol = ncol(I_theta),
       byrow = FALSE)
      #for the each kth row of single_pert_theta add to the kth row / kth column of I_theta
      I_theta <- I_theta - spt_wide - t(spt_wide)
      SE_CONVERGED <- TRUE
    }

    for (c in 1:ncol(I_theta))
    {
      pert_theta <- new_theta
      pert_theta[c] <- pert_theta[c] + h_N
      double_pert_theta <- sapply(X = seq(c, ncol(I_theta)),
        FUN = pl_theta,
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
      if(c < ncol(I_theta))
      {
        ## And fill the others in on the cth row/ column
        dpt[c, -(1:c)] <- dpt[-(1:c), c] <- double_pert_theta[-1]
      }

      I_theta <- I_theta + dpt
      # print(I_theta)
    }

    I_theta <- h_N ^ (- 2) * I_theta

    cov_theta <- tryCatch(expr = - solve(I_theta),
      error = function(err)
      {
        matrix(NA, nrow = nrow(I_theta), ncol = ncol(I_theta))
      }

      )
    # ------------------------- Estimate Cov(theta) using profile likelihood
    # if(any(diag(cov_theta) < 0)) {
    #   warning("Negative variance estimate. Increase the h_N_scale parameter and repeat variance estimation.")
    #   SE_CONVERGED <- FALSE
    # }

    se_theta <- tryCatch(expr = sqrt(diag(cov_theta)),
      warning = function(w)
      {
        matrix(NA, nrow = nrow(prev_theta))
        })
    if (any(is.na(se_theta)))
    {
      SE_CONVERGED <- FALSE
    }
    else
    {
      TRUE
    }

    return(list(coeff = data.frame(coeff = new_theta, se = se_theta),
                outcome_err_coeff = data.frame(coeff = new_gamma, se = NA),
                Bspline_coeff = cbind(k = comp_dat_val[, "k"], new_p),
                vcov = cov_theta,
                converged = CONVERGED,
                se_converged = SE_CONVERGED,
                converged_msg = CONVERGED_MSG,
                iterations = it,
                od_loglik_at_conv = od_loglik_theta))
  }
}

