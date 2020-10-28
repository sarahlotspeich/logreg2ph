hessian_row <- function(x, pm) {
  return(colSums(x * pm))
}
