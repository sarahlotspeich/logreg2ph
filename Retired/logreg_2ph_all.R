logreg_2ph <- function(Y_unval = NULL, Y_val = NULL, X_unval = NULL, X_val = NULL, C = NULL, interact = NULL,
         Validated = NULL, Bspline = NULL, data, initial_lr_params = "Zero",
         h_n_scale = NULL, noSE=FALSE, TOL = 1E-4, MAX_ITER = 1000, rescale = FALSE) {
  # Save sample sizes ---------------------------------------------------
  ## Phase - I ----------------------------------------------------------
  N <- nrow(data)
  ## Phase - II ---------------------------------------------------------
  n <- sum(data[, Validated])
  
  # Determine measurement error setting ---------------------------------
  ## Covariate error ----------------------------------------------------
  X_errors <- !is.null(X_unval)
  ## Outcome error ------------------------------------------------------
  Y_errors <- !is.null(Y_unval) 
  
  if (is.null(h_n_scale)) {
    if (rescale) {
      h_n_scale <- 4
    } else {
      h_n_scale <- 1
    }
  }
  
  # Reorder so that the n validated subjects are first ------------------
  data <- data[order(as.numeric(data[, Validated]), decreasing = TRUE), ]
  
  # Create interaction terms --------------------------------------------
  X_val_interact <- NULL
  if (length(interact) > 0) {
    for (int in 1:length(interact)) {
      interact_term <- data[, interact[[int]][1]]
      for (t in 2:length(interact[[int]])) {
        interact_term <- interact_term * data[, interact[[int]][t]]
      }
      data <- cbind(data, interact_term)
      colnames(data)[ncol(data)] <- paste0(interact[[int]][], collapse = "x")
      X_val_interact <- c(X_val_interact, paste0(interact[[int]][], collapse = "x"))
    }
  }
  
  # Standardize X_val, X_unval and C to N(0,1) ---------------------------
  re_data <- data
  if(rescale) {
    ## Shift by the sample mean ------------------------------------------
    re_shift <- c(0, as.numeric(colMeans(data[, c(X_val, X_val_interact, X_unval, C)], na.rm = TRUE)))
    ## Scale inversely by sample standard deviation ----------------------
    re_scale <- c(1, as.numeric(apply(data[, c(X_val, X_val_interact, X_unval, C)], MARGIN = 2, FUN = sd, na.rm = TRUE)))
    for (p in 1:length(c(X_val, X_val_interact, X_unval, C)))
    {
      re_data[, c(X_val, X_val_interact, X_unval, C)[p]] <- (re_data[,c(X_val, X_val_interact, X_unval, C)[p]] - re_shift[p+1]) / re_scale[p+1]
    }
  }
  
  if (Y_errors & !X_errors) {
    logreg_2ph_Yonly(N = N, n = n, Y_val = Y_val, Y_unval = Y_unval, 
                     X_val = X_val, C = C, X_val_interact = X_val_interact,
                     data = re_data, initial_lr_params = initial_lr_params, 
                     h_n_scale = h_n_scale, noSE = noSE, TOL = TOL, 
                     MAX_ITER = MAX_ITER, rescale = rescale)
  }
}