# 95\% Bootstrapped Interval for 95\% CI coverage of the SMLE 
boot_cover <- function(B = 1000, data, tru) {
  N <- nrow(data)
  cp <- vector(length = B)
  for (b in 1:B) {
    re_data <- data[sample(x = 1:N, size = N, replace = TRUE), ]
    cp[b] <- with(re_data, mean(beta1_est_smle - 1.96 * beta1_se_smle <= tru & tru <= beta1_est_smle + 1.96 * beta1_se_smle, na.rm = TRUE))
  }
  return(data.frame(mce_cp = sqrt(mean((cp - mean(cp))^2)), avg_cp = mean(cp), lb = quantile(cp, 0.025), ub = quantile(cp, 0.975)))
}

# 95\% Bootstrapped Interval for Bias of the SMLE 
boot_bias <- function(B = 1000, data, tru) {
  N <- nrow(data)
  bias <- vector(length = B)
  for (b in 1:B) {
    re_data <- data[sample(x = 1:N, size = N, replace = TRUE), ]
    bias[b] <- with(re_data, mean(beta1_est_smle - tru, na.rm = TRUE))
  }
  return(data.frame(mce_bias = sqrt(mean((bias - mean(bias))^2)), avg_bias = mean(bias), lb = quantile(bias, 0.025), ub = quantile(bias, 0.975)))
}