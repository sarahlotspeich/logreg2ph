library(magrittr)

sim_data <- read.csv("~sarahlotspeich/Downloads/VCCC_SIMULATED_DAT.csv", header=TRUE, stringsAsFactors=FALSE)

Y_val <- "ADE"; Y_unval <- "ADEstar"
contX_val <- "CD4"; contX_unval <- "CD4star"
binX_val <- c("PRIOR", "NAIVE"); binX_unval <- c("PRIORstar", "NAIVEstar")
Validated <- "V"

d <- length(c(contX_val, binX_val))
theta0 <- theta <- matrix(data = 0, nrow = (d+1), ncol = 1)
gamma0 <- gamma <- matrix(data = 0, nrow = (2*d + 1 + 1), ncol = 1)
lambda10 <- lambda1 <- matrix(data = 0, nrow = (d + 2 + 1), ncol = 1)
sigma0 <- sigma <- 1
lambda20 <- lambda2 <- matrix(data = 0, nrow = (d + 1 + 1), ncol = 1)
lambda30 <- lambda3 <- matrix(data = 0, nrow = (d + 1), ncol = 1)

source("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/TwoPhase_LogReg/logreg_2ph_MLE.R")

sim_data_complete <- setup_complete_data(dat = sim_data, 
                                         Y_unval="ADEstar", Y_val="ADE", 
                                         contX_unval="CD4star", contX_val="CD4", 
                                         binX_unval=c("PRIORstar", "NAIVEstar"), binX_val=c("PRIOR", "NAIVE"), 
                                         Validated = "V_CC")

val_id <- which(sim_data[,"V_CC"])
n <- sum(sim_data[,"V_CC"])
m <- nrow(unique(sim_data[val_id,c("CD4", "PRIOR", "NAIVE")]))

system.time(expr = cc_mle_res <- nlm(
  f = observed_data_loglik_fp, p = c(theta, gamma, lambda1, lambda2,lambda3, sigma), hessian = TRUE, 
  n = nrow(sim_data), n_v = n, m = m, Y_unval="ADEstar", Y_val="ADE", 
  contX_unval="CD4star", contX_val="CD4", binX_unval=c("PRIORstar", "NAIVEstar"), binX_val=c("PRIOR", "NAIVE"), 
  C=NULL, comp_dat_all = sim_data_complete, 
  print.level = 2, iterlim = 1000
))

nr_mle <- function()
{
  alpha <- 0
  beta <- 0
  Delta_l <- Inf
  l <- log_likelihood(x, y, beta, alpha)
  iter <- 0
  while(abs(Delta_l) > TOL & iter < MAX_ITER)
  {
    iter <- iter + 1
    g <- gradient(x, y, beta, alpha)
    hess <- hessian(x, y, beta, alpha)
    H_inv <- solve(hess)
    
    Delta <- H_inv %*% g
    Delta_beta <- Delta[1]
    Delta_alpha <- Delta[2]
    
    alpha <- alpha + Delta_alpha
    beta <- beta + Delta_beta
    
    l_new <- log_likelihood(x, y, beta, alpha)
    Delta_l <- l - l_new
    l <- l_new
  }
  return(matrix(c(alpha, beta), nrow =  2))
}

gradient_logreg_2ph_MLE(theta = theta0, gamma = gamma0, lambda1 = lambda10, lambda2 = lambda20, lambda3 = lambda30, sigma = sigma0, 
                        n = n, m = m, complete_data = sim_data_complete, 
                        Y_unval="ADEstar", Y_val="ADE", 
                        contX_unval="CD4star", contX_val="CD4", 
                        binX_unval=c("PRIORstar", "NAIVEstar"), binX_val=c("PRIOR", "NAIVE"))

cc_mle_res <- nlm(
  f = observed_data_loglik_fp, p = c(theta, gamma, lambda1, lambda2,lambda3, sigma), hessian = TRUE, 
  n = nrow(sim_data), n_v = n, m = m, Y_unval="ADEstar", Y_val="ADE", 
  contX_unval="CD4star", contX_val="CD4", binX_unval=c("PRIORstar", "NAIVEstar"), binX_val=c("PRIOR", "NAIVE"), 
  C=NULL, comp_dat_all = sim_data_complete, 
  print.level = 2, iterlim = 1
)
