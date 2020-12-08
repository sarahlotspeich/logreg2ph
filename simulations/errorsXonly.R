# If running on laptop:
library(logreg2ph)
#source("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/Paper - SMLE logistic regression/Tang_twophase_loglik_binaryX.R")

# If running on ACCRE:
# libFolder = "~/R/rlib-3.4.3"
# library(logreg2ph, lib.loc = libFolder)
# source("/home/lotspesc/SMLE_LogReg/Tang_twophase_loglik_binaryX.R")

audit <- "Validated case-control"

sim_seed <- 918
num_sims <- 1000

# True parameter values for P(Y|X,Z) ----------------------
beta0 <- -1; beta1 <- 1; beta2 <- -0.5

N <- 1000
pv <- 0.1
n <- N * pv
set.seed(sim_seed)

smle_res <- data.frame(sim = 1:num_sims, N, n, beta0, beta1, beta2, PhaseII = audit,
                       beta0_est = NA, beta1_est = NA, beta2_est = NA,
                       beta0_se = NA, beta1_se = NA, beta2_se = NA,
                       conv = NA, se_conv = NA, conv_msg = NA)

for(s in 1:num_sims)
{
  # Generate true values Y, X, Z ----------------------------
  Z <- rbinom(n = N, size = 1, prob = 0.25)
  X <- rnorm(n = N, mean = 0, sd = 1)
  Y <- rbinom(n = N, size = 1,
              prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z))) ^ (- 1))
  
  # Generate error-prone X* = X + U -------------------------
  U <- rnorm(n = N, mean = 0, sd = sqrt(0.1))
  Xstar <- X + U
  
  # Audit
  if(audit == "SRS")
  {
    V <- seq(1, N) %in% sample(x = seq(1, N), size = n, replace = FALSE)
  }
  if(audit == "Validated case-control")
  {
    V <- seq(1, N) %in%
      c(sample(x = which(Y == 0), size = 0.5 * n, replace = FALSE),
        sample(x = which(Y == 1), size = 0.5 * n, replace = FALSE))
  }
  
  ### Construct B-spline basis -------------------------------
  ### We chose cubic B-splines, with 20 df for N = 1000 ------
  nsieve <- 20
  B <- matrix(0, nrow = N, ncol = nsieve)
  B[which(Z == 0),1:(0.75 * nsieve)] <- splines::bs(x = Xstar[which(Z == 0)], df = 0.75 * nsieve, Boundary.knots = range(Xstar[which(Z == 0)]), intercept = TRUE)
  B[which(Z == 1),(0.75 * nsieve + 1):nsieve] <- splines::bs(x = Xstar[which(Z == 1)], df = 0.25 * nsieve, Boundary.knots = range(Xstar[which(Z == 1)]), intercept = TRUE)
  colnames(B) <- paste0("bs", seq(1, nsieve))
  
  # Build dataset --------------------------------------------
  sdat <- cbind(Y, X, Xstar, Z, V, B)
  # Make Phase-II variables Y, X NA for unaudited subjects ---
  sdat[!V, "X"] <- NA
  
  suppressMessages(
    smle <- logreg2ph::logreg2ph_rw(Y_unval = NULL,
                                    Y_val = "Y",
                                    X_unval = "Xstar",
                                    X_val = "X",
                                    C = "Z",
                                    Validated = "V",
                                    Bspline = colnames(B),
                                    data = sdat,
                                    noSE = FALSE,
                                    MAX_ITER = 1000,
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
  
  if (s %%10 == 0) {print(paste(s, "simulations complete."))}
  
  write.csv(smle_res,
            paste0("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/logreg_2ph/test_errorsXonly_N", N, "_n", n, ".csv"),
            row.names = FALSE)
}
