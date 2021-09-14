nrst <- function(target, Bspline_coeff, k = 2) {
  target_rep <- do.call(rbind, replicate(n = nrow(Bspline_coeff), expr = target, simplify = FALSE))
  neighborhood <- Bspline_coeff[, 1:(ncol(Bspline_coeff) / 2)]
  coeff <- Bspline_coeff[, -c(1:(ncol(Bspline_coeff) / 2))]
  diff <- rowSums(abs(target_rep - neighborhood))
  coeff_ord <- coeff[as.vector(order(diff)), ]
  return(colMeans(coeff_ord[1:k, ]))
}
