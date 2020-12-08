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
theta0 <- -log(specY / (1 - specY))
theta1 <- -theta0 - log((1 - sensY) / sensY)
theta2 <- -0.2
theta3 <- -0.2
theta4 <- -0.1

sensX <- 0.75
specX <- 0.75
delta0 <- -log(specX/(1-specX))
delta1 <- -delta0 - log((1-sensX)/sensX)
delta2 <- -0.2

gamma0 <- 0

N <- 1000
pv <- 0.25
n <- N * pv
set.seed(sim_seed)
deg <- 1

for (nsieve in c(18, 20, 22, 24, 26)) {
  smle_res <- data.frame(sim = 1:num_sims, N, n, beta0, beta1, beta2, theta5, PhaseII = audit,
                         deg, nsieve, beta0_est = NA, beta1_est = NA, beta2_est = NA,
                         beta0_se = NA, beta1_se = NA, beta2_se = NA,
                         conv = NA, se_conv = NA, conv_msg = NA)

  for(s in 1:num_sims)
  {
    # Generate Z ~ Norm(0, 1)
    Z <- rnorm(n = N, mean = 0, sd = 1)
    X <- rbinom(n = N, size = 1, prob = (1 + exp(-(gamma0)))^(-1))
    Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z)))^(-1))
    Xstar <- rbinom(n = N, size = 1, prob = (1 + exp(-(delta0 + delta1 * X + delta2 * Z)))^(-1))
    Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(-(theta0 + theta1 * Y + theta2 * X + theta3 * Xstar + theta4 * Z + theta5 * Z ^ 2)))^(-1))

    # Audit
    if(audit == "SRS")
    {
      V <- seq(1,N) %in% sample(x = seq(1,N), size = n, replace = FALSE)
    }
    if(audit == "Unvalidated case-control")
    {
      V <- seq(1,N) %in% c(sample(x = which(Ystar == 0), size = 0.5*n, replace = FALSE), sample(x = which(Ystar == 1), size = 0.5*n, replace = FALSE))
    }

    B <- matrix(0, nrow = N, ncol = nsieve)
    B[which(Xstar == 0), 1:(nsieve / 2)] <- splines::bs(x = Z[which(Xstar == 0)], df = nsieve / 2, degree = deg, Boundary.knots = range(Z[which(Xstar == 0)]), intercept = TRUE)
    B[which(Xstar == 1), (nsieve / 2 + 1):nsieve] <- splines::bs(x = Z[which(Xstar == 1)], df = nsieve / 2, degree = deg, Boundary.knots = range(Z[which(Xstar == 1)]), intercept = TRUE)
    colnames(B) <- paste0("bs", seq(1, nsieve))

    sdat <- cbind(Y, X, Ystar, Xstar, Z, Z2 = Z ^ 2, V, B)
    sdat[!V, c("Y", "X")] <- NA

    suppressMessages(
      smle <- logreg2ph::logreg2ph_rw(Y_unval = "Ystar", Y_val = "Y", X_unval = "Xstar", X_val = "X", C = "Z",
                                      Validated = "V", Bspline = colnames(B), data = sdat,
                                      gamma_pred = c("Xstar", "Y", "X", "Z", "Z2"),
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

    if (s %%10 == 0) {print(paste(s, "simulations complete."))}

    write.csv(smle_res, paste0("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/logreg_2ph/tuneBspline/misspecYstarmod_correctSMLE_deg", deg, "_", nsieve, "sieves_seed", sim_seed, "_", ifelse(audit == "SRS", "SRS", "CC"), ".csv"), row.names = FALSE)
  }
}
