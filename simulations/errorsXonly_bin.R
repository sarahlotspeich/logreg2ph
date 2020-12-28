source('~/Dropbox/Vanderbilt/Lotspeich_Dissertation/logreg_2ph/R/logreg2ph_Xonly.R')
source('~/Dropbox/Vanderbilt/Lotspeich_Dissertation/logreg_2ph/R/logreg2ph_Xonly_nlm.R')
source('~/Dropbox/Vanderbilt/Lotspeich_Dissertation/logreg_2ph/R/observed_data_loglik_Xonly.R')

audit <- "SRS"

sim_seed <- 918
prevZ <- 0.25
num_sims <- 1000

beta0_Xb <- -0.65
beta1_Xb <- -0.20
beta2_Xb <- -0.10

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
delta2 <- 0 #-0.2
delta3 <- 0.5

gamma0 <- 0

N <- 1000
#pv <- 0.1

smle_res <- data.frame(sim = 1:num_sims,
                       beta0_em = NA, beta1_em = NA, beta2_em = NA,
                       beta0_nlm = NA, beta1_nlm = NA, beta2_nlm = NA,
                       conv_em = NA, conv_nlm = NA,
                       min_od_em = NA, min_od_nlm = NA)
for (pv in c(0.1, 0.25, 0.5, 0.75)) {
  n <- N * pv
  set.seed(sim_seed)
  for(s in 1:num_sims)
  {
    # Generate true values Y, X, Z ----------------------------
    Z <- rbinom(n = N, size = 1, prob = prevZ)
    X <- rbinom(n = N, size = 1, prob = (1 + exp(-(gamma0)))^(-1))
    Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(beta0_Xb + beta1_Xb*X + beta2_Xb*Z)))^(-1))

    # Generate error-prone X* = X + U -------------------------
    Xstar <- rbinom(n = N, size = 1, prob = (1 + exp(-(delta0+delta1*X+delta2*Y+delta3*Z)))^(-1))

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
    nsieve <- 4
    B <- matrix(0, nrow = N, ncol = nsieve)
    B[which(Z == 0 & Xstar == 0),1] <- 1
    B[which(Z == 0 & Xstar == 1),2] <- 1
    B[which(Z == 1 & Xstar == 0),3] <- 1
    B[which(Z == 1 & Xstar == 1),4] <- 1
    colnames(B) <- paste0("bs", seq(1, nsieve))

    # Build dataset --------------------------------------------
    sdat <- cbind(Y, X, Xstar, Z, V, B)
    # Make Phase-II variables Y, X NA for unaudited subjects ---
    sdat[!V, "X"] <- NA
    suppressMessages(
      em <- logreg2ph_Xonly(Y_val = "Y",
                            X_unval = "Xstar",
                            X_val = "X",
                            C = "Z",
                            Validated = "V",
                            Bspline = colnames(B),
                            data = sdat,
                            theta_pred = NULL,
                            initial_lr_params = "Zero",
                            h_N_scale = 1,
                            noSE = FALSE,
                            TOL = 1E-4,
                            MAX_ITER = 1000)
    )

    smle_res$beta0_em[s] <- em$Coefficients$Coefficient[1]
    smle_res$beta1_em[s] <- em$Coefficients$Coefficient[2]
    smle_res$beta2_em[s] <- em$Coefficients$Coefficient[3]
    smle_res$conv_em[s] <- em$converged
    smle_res$min_od_em[s] <- em$od_loglik_at_conv

    # Check that the parameters are actually maximizing
    # em$all_iter_od_loglik %>%
    #   data.frame() %>% ggplot(aes(x = iter, y = od_ll)) + geom_line()
    suppressMessages(
      nlm <- logreg2ph_Xonly_nlm(Y_val = "Y",
                                 X_unval = "Xstar",
                                 X_val = "X",
                                 C = "Z",
                                 Validated = "V",
                                 Bspline = colnames(B),
                                 data = sdat,
                                 theta_pred = NULL,
                                 initial_lr_params = "Zero",
                                 TOL = 1E-4,
                                 MAX_ITER = 1000)
    )
    smle_res$beta0_nlm[s] <- nlm$Coefficients$Coefficient[1]
    smle_res$beta1_nlm[s] <- nlm$Coefficients$Coefficient[2]
    smle_res$beta2_nlm[s] <- nlm$Coefficients$Coefficient[3]
    smle_res$conv_nlm[s] <- nlm$converged
    smle_res$min_od_nlm[s] <- nlm$od_loglik_at_conv

    if (s %%10 == 0) {print(paste(s, "simulations complete."))}

    write.csv(smle_res,
              paste0("~sarahlotspeich/Downloads/test_errorsXonly_N", N, "_n", n, ".csv"),
              row.names = FALSE)
  }
}

d <- read.csv("~sarahlotspeich/Downloads/test_errorsXonly_N1000_n100.csv") %>%
  dplyr::mutate(n = 100) %>%
  dplyr::bind_rows(read.csv("~sarahlotspeich/Downloads/test_errorsXonly_N1000_n250.csv") %>%
                     dplyr::mutate(n = 250)) %>%
  dplyr::bind_rows(read.csv("~sarahlotspeich/Downloads/test_errorsXonly_N1000_n500.csv") %>%
                     dplyr::mutate(n = 500)) %>%
  dplyr::bind_rows(read.csv("~sarahlotspeich/Downloads/test_errorsXonly_N1000_n750.csv") %>%
                     dplyr::mutate(n = 750))

d %>% dplyr::group_by(n) %>%
  dplyr::select(-sim) %>%
  dplyr::summarize_all(mean, na.rm = TRUE)
