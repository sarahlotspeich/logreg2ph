# If running on laptop:
library(logreg2ph)
#source("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/Paper - SMLE logistic regression/Tang_twophase_loglik_binaryX.R")

# If running on ACCRE:
# libFolder = "~/R/rlib-3.4.3"
# library(logreg2ph, lib.loc = libFolder)
# source("/home/lotspesc/SMLE_LogReg/Tang_twophase_loglik_binaryX.R")

audit <- "Unvalidated case-control"

sim_seed <- 918

num_sims <- 1000

# True parameter values for P(Y|X,Z) ----------------------
beta0 <- -1; beta1 <- 1; beta2 <- -0.5

# Parameters for error model P(Y*|X*,Y,X,Z) ---------------
sensY <- 0.95; specY <- 0.90
theta0 <- - log(specY / (1 - specY))
theta1 <- - theta0 - log((1 - sensY) / sensY)
theta2 <- 1; theta3 <- - 0.5

for (N in c(1000, 2000)) {
  for (pv in c(0.1, 0.25, 0.4)) {
    n <- N * pv
    set.seed(sim_seed)

    smle_res <- data.frame(sim = 1:num_sims, N, n, beta0, beta1, beta2, PhaseII = audit,
                           beta0_est = NA, beta1_est = NA, beta2_est = NA,
                           beta0_se = NA, beta1_se = NA, beta2_se = NA,
                           conv = NA, se_conv = NA, conv_msg = NA, iter = NA)

    for(s in 1:num_sims)
    {
      # Generate true values Y, X, Z ----------------------------
      Z <- rbinom(n = N, size = 1, prob = 0.25)
      X <- rnorm(n = N, mean = 0, sd = 1)
      Y <- rbinom(n = N, size = 1,
                  prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z))) ^ (- 1))

      # Generate error-prone Y* from error model P(Y*|Y,X,Z) --
      Ystar <- rbinom(n = N, size = 1,
                      prob = (1 + exp(- (theta0 + theta1 * Y + theta2 * X + theta3 * Z))) ^ (- 1))

      # Audit
      if(audit == "SRS")
      {
        V <- seq(1, N) %in% sample(x = seq(1, N), size = n, replace = FALSE)
      }
      if(audit == "Unvalidated case-control")
      {
        V <- seq(1, N) %in%
          c(sample(x = which(Ystar == 0), size = 0.5 * n, replace = FALSE),
            sample(x = which(Ystar == 1), size = 0.5 * n, replace = FALSE))
      }

      sdat <- cbind(Y, X, Ystar, Z, V)
      sdat[!V, c("Y")] <- NA

      suppressMessages(
        smle <- logreg2ph::logreg2ph_rw(Y_unval = "Ystar",
                                        Y_val = "Y",
                                        C = c("X", "Z"),
                                        Validated = "V",
                                        data = sdat,
                                        noSE = FALSE,
                                        MAX_ITER = 10000,
                                        TOL = 1E-4)
      )

      smle_res$beta0_est[s] <- smle$Coefficients$Coefficient[1]
      smle_res$beta1_est[s] <- smle$Coefficients$Coefficient[2]
      smle_res$beta2_est[s] <- smle$Coefficients$Coefficient[3]
      smle_res$beta0_se[s] <- smle$Coefficients$SE[1]
      smle_res$beta1_se[s] <- smle$Coefficients$SE[2]
      smle_res$beta2_se[s] <- smle$Coefficients$SE[3]
      smle_res$conv[s] <- smle$converged
      smle_res$conv_msg[s] <- smle$converged_msg
      smle_res$iter[s] <- smle$iterations

      if (s %%10 == 0) {print(paste(s, "simulations complete."))}

      write.csv(smle_res,
                paste0("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/logreg_2ph/errorsY_only_res/test_errorsYonly_N", N, "_n", n, ".csv"),
                row.names = FALSE)
    }
  }
}
