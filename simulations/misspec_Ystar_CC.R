# If running on laptop:
library(logreg2ph)
source("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/Paper - SMLE logistic regression/Tang_twophase_loglik_binaryX.R")

# If running on ACCRE:
# libFolder = "~/R/rlib-3.4.3"
# library(logreg2ph, lib.loc = libFolder)
# source("/home/lotspesc/SMLE_LogReg/Tang_twophase_loglik_binaryX.R")

audit <- "Unvalidated case-control"
theta5 <- 0

args <- commandArgs(TRUE)

sim_seed <- 918 #as.integer(args)

num_sims <- 1000

beta0 <- -0.65
beta1 <- -0.20
beta2 <- -0.10

sensY <- 0.95
specY <- 0.90
theta0 <- -log(specY/(1-specY))
theta1 <- -theta0 - log((1-sensY)/sensY)
theta2 <- -0.2
theta3 <- -0.2
theta4 <- -0.1

sensX <- 0.75
specX <- 0.75
delta0 <- -log(specX/(1-specX))
delta1 <- -delta0 - log((1-sensX)/sensX)
delta2 <- -0.2
delta3 <- -0.1

gamma0 <- 0

N <- 1000
pv <- 0.25
n <- N * pv
set.seed(sim_seed)

sim_res <- data.frame(sim = 1:num_sims, N, n, beta0, beta1, beta2, theta5, PhaseII = audit,
                      prevY = NA, prevYstar = NA, mcY = NA, fprY = NA, fnrY = NA, mcX = NA, fprX = NA, fnrX = NA,
                      beta0_est = NA, beta1_est = NA, beta2_est = NA, beta0_se = NA, beta1_se = NA, beta2_se = NA,
                      conv = NA, se_conv = NA, conv_msg = NA)

naive_res <- complete_res <- tang_res <- ht_res <- smle_res <- sim_res
naive_res$method <- "Naive"
complete_res$method <- "Complete-data"
tang_res$method <- "MLE"
smle_res$method <- "SMLE"
ht_res$method <- "HT"

for(s in 1:num_sims)
{
  # Generate Z ~ Norm(0, 1)
  Z <- rnorm(n = N, mean = 0, sd = 1)
  X <- rbinom(n = N, size = 1, prob = (1 + exp(-(gamma0)))^(-1))
  Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z)))^(-1))
  Xstar <- rbinom(n = N, size = 1, prob = (1 + exp(-(delta0 + delta1 * X + delta2 * Z + delta3 * Z ^ 2)))^(-1))
  Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(-(theta0 + theta1 * Y + theta2 * X + theta3 * Xstar + theta4 * Z + theta5 * Z ^ 2)))^(-1))

  smle_res$prevY[s] <- mean(Y)
  smle_res$prevYstar[s] <- mean(Ystar)

  smle_res$mcY[s] <- mean(Y != Ystar)
  smle_res$fprY[s] <- mean(Ystar[Y == 0])
  smle_res$fnrY[s] <- 1-mean(Ystar[Y==1])
  smle_res$mcX[s] <- mean(X != Xstar)
  smle_res$fprX[s] <- mean(Xstar[X == 0])
  smle_res$fnrX[s] <- 1-mean(Xstar[X == 1])

  # naive model
  naive <- glm(Ystar ~ Xstar + Z, family = "binomial")
  naive_res$beta0_est[s] <- naive$coefficients[1]
  naive_res$beta1_est[s] <- naive$coefficients[2]
  naive_res$beta2_est[s] <- naive$coefficients[3]
  naive_res$beta0_se[s] <- sqrt(diag(vcov(naive)))[1]
  naive_res$beta1_se[s] <- sqrt(diag(vcov(naive)))[2]
  naive_res$beta2_se[s] <- sqrt(diag(vcov(naive)))[3]

  # Audit
  if(audit == "SRS")
  {
    V <- seq(1,N) %in% sample(x = seq(1,N), size = n, replace = FALSE)
  }
  if(audit == "Unvalidated case-control")
  {
    V <- seq(1,N) %in% c(sample(x = which(Ystar == 0), size = 0.5*n, replace = FALSE), sample(x = which(Ystar == 1), size = 0.5*n, replace = FALSE))
  }

  cc <- glm(Y[V] ~ X[V] + Z[V], family = "binomial")
  complete_res$beta0_est[s] <- cc$coefficients[1]
  complete_res$beta1_est[s] <- cc$coefficients[2]
  complete_res$beta2_est[s] <- cc$coefficients[3]
  complete_res$beta0_se[s] <- sqrt(diag(vcov(cc)))[1]
  complete_res$beta1_se[s] <- sqrt(diag(vcov(cc)))[2]
  complete_res$beta2_se[s] <- sqrt(diag(vcov(cc)))[3]

  # IPW complete-data model/Horvitz Thompson estimator
  ht <- glm(Y[V] ~ X[V] + Z[V], family = "binomial",
            weights = ifelse(Ystar[V] == 0,
                             1/((0.5*n)/(table(Ystar)[1])),
                             1/((0.5*n)/(table(Ystar)[2]))))
  ht_res$beta0_est[s] <- ht$coefficients[1]
  ht_res$beta1_est[s] <- ht$coefficients[2]
  ht_res$beta2_est[s] <- ht$coefficients[3]
  ht_res$beta0_se[s] <- sqrt(diag(vcov(ht)))[1]
  ht_res$beta1_se[s] <- sqrt(diag(vcov(ht)))[2]
  ht_res$beta2_se[s] <- sqrt(diag(vcov(ht)))[3]

  nsieve <- 24
  B <- matrix(0, nrow = N, ncol = nsieve)
  B[which(Xstar == 0), 1:12] <- splines::bs(x = Z[which(Xstar == 0)], df = 12, Boundary.knots = range(Z[which(Xstar == 0)]), intercept = TRUE)
  B[which(Xstar == 1), 13:24] <- splines::bs(x = Z[which(Xstar == 1)], df = 12, Boundary.knots = range(Z[which(Xstar == 1)]), intercept = TRUE)
  colnames(B) <- paste0("bs", seq(1, nsieve))

  sdat <- cbind(Y, X, Ystar, Xstar, Z, V, B)
  sdat[!V,c("Y", "X")] <- NA

  suppressMessages(
    smle <- logreg2ph::logreg2ph_rw(Y_unval = "Ystar", Y_val = "Y", X_unval = "Xstar", X_val = "X", C = "Z",
                                    Validated = "V", Bspline = colnames(B), data = sdat,
                                    noSE = FALSE, MAX_ITER = 1000, TOL = 1E-4)
  )

  smle_res$beta0_est[s] <- smle$Coefficients$Coefficient[1]
  smle_res$beta1_est[s] <- smle$Coefficients$Coefficient[2]
  smle_res$beta2_est[s] <- smle$Coefficients$Coefficient[3]
  smle_res$beta0_se[s] <- smle$Coefficients$SE[1]
  smle_res$beta1_se[s] <- smle$Coefficients$SE[2]
  smle_res$beta2_se[s] <- smle$Coefficients$SE[3]
  smle_res$conv[s] <- smle$converged
  smle_res$conv_msg[s] <- smle$converged_msg

  # Fully-parametric likelihood (Tang et al. 2015) -----
  suppressWarnings(fit_Tang <- nlm(Tang_twophase_loglik, params <- rep(0,12), hessian=TRUE, Val = "V", Y_unval="Ystar", Y_val="Y", X_unval="Xstar", X_val="X", C = "Z", data = sdat))
  tang_res$beta0_est[s] <- ifelse(fit_Tang$code == 1, fit_Tang$estimate[9], NA)
  tang_res$beta0_se[s] <- tryCatch(expr = sqrt(diag(solve(fit_Tang$hessian)))[9],  error = function(err) {NA})
  tang_res$beta1_est[s] <- ifelse(fit_Tang$code == 1, fit_Tang$estimate[10], NA)
  tang_res$beta1_se[s] <- tryCatch(expr = sqrt(diag(solve(fit_Tang$hessian)))[10],  error = function(err) {NA})
  tang_res$beta2_est[s] <- ifelse(fit_Tang$code == 1, fit_Tang$estimate[11], NA)
  tang_res$beta2_se[s] <- tryCatch(expr = sqrt(diag(solve(fit_Tang$hessian)))[11], error = function(err) {NA})
  tang_res$conv[s] <- fit_Tang$code <= 2
  # ----- Fully-parametric likelihood (Tang et al. 2015)

  if (s %%10 == 0) {print(paste(s, "simulations complete."))}

  write.csv(rbind(smle_res, tang_res, naive_res), paste0("misspecYstarmod_theta", theta5, "_", N, "_pv", 100*pv,"_seed", sim_seed, "_", ifelse(audit == "SRS", "SRS", "CC"), ".csv"), row.names = FALSE)
}
