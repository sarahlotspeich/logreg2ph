#' Observed-data log-likelihood
#'
#' This function returns the value of the observed-data log-likelihood (equation (2) in Lotspeich et al. (2020))
#' for a given dataset and parameter values `theta`, `gamma`, and `p`.
#
#'
#' @param N Phase I sample size
#' @param n Phase II sample size
#' @param Y_unval Column with the unvalidated outcome (can be name or numeric index)
#' @param Y_val Column with the validated outcome (can be name or numeric index)
#' @param X_unval Column(s) with the unvalidated predictors (can be name or numeric index)
#' @param X_val Column(s) with the validated predictors (can be name or numeric index)
#' @param C (Optional) Column(s) with additional error-free covariates (can be name or numeric index)
#' @param Bspline Vector of columns containing the B-spline basis functions (can be name or numeric index)
#' @param comp_dat_all Augmented dataset containing rows for each combination of unvalidated subjects' data with values from Phase II (a matrix)
#' @param theta Parameters for the analysis model (a column vector)
#' @param gamma Parameters for the outcome error model (a column vector)
#' @param p B-spline coefficients for the approximated covariate error model (a matrix)
#' @return Scalar value of the function

observed_data_loglik <- function(N, n, Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL,
                                 Bspline = NULL, comp_dat_all, theta, gamma, p) {
  sn <- ncol(p)
  m <- nrow(p)

  # For validated subjects --------------------------------------------------------
  #################################################################################
  ## Sum over log[P_theta(Yi|Xi)] -------------------------------------------------
  pY_X <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n), c(X_val, C)]) %*% theta)))
  pY_X <- ifelse(as.vector(comp_dat_all[c(1:n), c(Y_val)]) == 0, 1 - pY_X, pY_X)
  return_loglik <- sum(log(pY_X))
  ## ------------------------------------------------- Sum over log[P_theta(Yi|Xi)]
  #################################################################################
  ## Sum over log[P(Yi*|Xi*,Yi,Xi)] -----------------------------------------------
  pYstar <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[c(1:n), c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar <- ifelse(as.vector(comp_dat_all[c(1:n), Y_unval]) == 0, 1 - pYstar, pYstar)
  return_loglik <- return_loglik + sum(log(pYstar))
  ## ----------------------------------------------- Sum over log[P(Yi*|Xi*,Yi,Xi)]
  #################################################################################
  ## Sum over I(Xi=xk)Bj(Xi*)log p_kj ---------------------------------------------
  pX <- p[comp_dat_all[c(1:n), "k"], ]
  log_pX <- log(pX)
  log_pX[log_pX == -Inf] <- 0
  return_loglik <- return_loglik + sum(comp_dat_all[c(1:n), Bspline] * log_pX)
  ## --------------------------------------------- Sum over I(Xi=xk)Bj(Xi*)log q_kj
  #################################################################################
  # -------------------------------------------------------- For validated subjects

  # For unvalidated subjects ------------------------------------------------------
  ## Calculate P_theta(y|x) for all (y,xk) ----------------------------------------
  pY_X <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n), c(X_val, C)]) %*% theta)))
  pY_X[which(comp_dat_all[-c(1:n), Y_val] == 0)] <- 1-pY_X[which(comp_dat_all[-c(1:n), Y_val] == 0)]
  ## ---------------------------------------- Calculate P_theta(y|x) for all (y,xk)
  ################################################################################
  ## Calculate P(Yi*|Xi*,y,xk) for all (y,xk) ------------------------------------
  pYstar <- 1 / (1 + exp(-as.numeric(cbind(int = 1, comp_dat_all[-c(1:n), c(X_unval, Y_val, X_val, C)]) %*% gamma)))
  pYstar[which(comp_dat_all[-c(1:n), Y_unval] == 0)] <- 1 - pYstar[which(comp_dat_all[-c(1:n), Y_unval] == 0)]
  ## ------------------------------------ Calculate P(Yi*|Xi*,y,xk) for all (y,xk)
  ################################################################################
  ## Calculate Bj(Xi*) p_kj for all (k,j) ----------------------------------------
  pX <- p[comp_dat_all[-c(1:n), "k"], ]
  ## ---------------------------------------- Calculate Bj(Xi*) p_kj for all (k,j)
  ################################################################################
  ## Calculate P(y|xk) x P(Y*|X*,y,xk) x Bj(X*) x p_kj ---------------------------
  person_sum <- rowsum(c(pY_X * pYstar) * comp_dat_all[-c(1:n), Bspline] * pX, group = rep(seq(1, (N - n)), times = 2 * m))
  person_sum <- rowSums(person_sum)
  log_person_sum <- log(person_sum)
  log_person_sum[log_person_sum == -Inf] <- 0
  ## And sum over them all -------------------------------------------------------
  return_loglik <- return_loglik + sum(log_person_sum)
  ################################################################################
  # ----------------------------------------------------- For unvalidated subjects
  return(return_loglik)
}
