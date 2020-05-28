source("~sarahlotspeich/Dropbox/Vanderbilt/Lotspeich_Dissertation/TwoPhase_LogReg/logreg_2ph_MLE.R")

library(magrittr); library(splines); library(dplyr)

VCCC_dat <- read.csv("~sarahlotspeich/Box/vccc-validation-data/VCCC_analysis_dat_V3.csv", header=TRUE, stringsAsFactors=FALSE) %>% 
  dplyr::mutate(SQRT_CD4_VAL = sqrt(CD4_COUNT_VAL), SQRT_CD4_UNVAL = sqrt(CD4_COUNT_UNVAL),
                SQRT_CD4_VAL10 = sqrt(CD4_COUNT_VAL/10), SQRT_CD4_UNVAL10 = sqrt(CD4_COUNT_UNVAL/10),
                PRIOR_NAIVE_UNVAL = paste(ifelse(PRIOR_ADE_UNVAL == 0, "No prior ADE", "Prior ADE"),ifelse(ART_NAIVE_ENROLL_UNVAL == 0, "Not ART naive", "ART naive"), sep = "/")
  )

Y_unval <- "ADE_WITHIN1YR_UNVAL"; Y_val <- "ADE_WITHIN1YR_VAL"
contX_unval <- "SQRT_CD4_UNVAL10"; contX_val <- "SQRT_CD4_VAL10"
binX_unval <- c("PRIOR_ADE_UNVAL", "ART_NAIVE_ENROLL_UNVAL"); binX_val <- c("PRIOR_ADE_VAL", "ART_NAIVE_ENROLL_VAL")
C <- NULL
Validated <- "V"

vccc_mle_sim <- function(sim_seed, controls_mult, num_sims, data)
{
  set.seed(sim_seed)
  for(nv in c((1+controls_mult)*with(VCCC_dat, table(ADE_WITHIN1YR_UNVAL))[2]))
  {
    sim_res <- data.frame(rep = 1:num_sims, Method = "MLE", PhaseII = "Case-control", nv, 
                          alpha = NA, alphaSE = NA, betaCD4 = NA, betaCD4SE = NA, 
                          betaPriorADE = NA, betaPriorADESE = NA, betaARTNaive = NA, betaARTNaiveSE = NA,
                          converged = NA, converged_message = NA, od_loglik_theta = NA,
                          audit_resampled = FALSE)
    
    for (r in 1:num_sims)
    {
      V <- seq(1, nrow(VCCC_dat)) %in% c(which(VCCC_dat$ADE_WITHIN1YR_UNVAL == 1), sample(which(VCCC_dat$ADE_WITHIN1YR_UNVAL == 0), size = controls_mult*with(VCCC_dat, table(ADE_WITHIN1YR_UNVAL))[2], replace = FALSE))
      
      nv_sample <- VCCC_dat
      nv_sample[!V, c("ADE_WITHIN1YR_VAL", "PRIOR_ADE_VAL", "SQRT_CD4_VAL10", "ART_NAIVE_ENROLL_VAL")] <- NA
      nv_sample %<>% cbind(V) 
      
      n <- nrow(nv_sample)
      n_v <- sum(nv_sample[,Validated])
      
      # Reorder so that the n_v validated subjects are first ------------
      nv_sample <- nv_sample[order(as.numeric(nv_sample[,Validated]), decreasing = TRUE),]
      
      # Save distinct X -------------------------------------------------
      x_obs <- data.frame(unique(nv_sample[1:n_v,c(contX_val, binX_val)]))
      x_obs <- data.frame(x_obs[order(x_obs[,1]),])
      colnames(x_obs) <- c(contX_val, binX_val)
      m <- nrow(x_obs)
      x_obs_stacked <- do.call(rbind, replicate(n = (n-n_v), expr = x_obs, simplify = FALSE))
      x_obs_stacked <- data.frame(x_obs_stacked[order(x_obs_stacked[,1]),])
      colnames(x_obs_stacked) <- c(contX_val, binX_val)
      
      # Save static (X*,Y*,X,Y,C) since they don't change ---------------
      comp_dat_val <- nv_sample[c(1:n_v),c(Y_unval, contX_unval, binX_unval, C, contX_val, binX_val, Y_val)]
      comp_dat_val <- merge(x = comp_dat_val, y = data.frame(x_obs, k = 1:m), all.x = TRUE)
      comp_dat_val <- comp_dat_val[,c(Y_unval, contX_unval, binX_unval, C, contX_val, binX_val, Y_val, "k")]
      comp_dat_val <- data.matrix(comp_dat_val)
      # 2 (m x n)xd matrices (y=0/y=1) of each (one column per person, --
      # one row per x) --------------------------------------------------
      suppressWarnings(comp_dat_unval <- cbind(nv_sample[-c(1:n_v),c(Y_unval, contX_unval, binX_unval, C)],
                                               x_obs_stacked))
      comp_dat_y0 <- data.frame(comp_dat_unval, Y = 0)
      comp_dat_y1 <- data.frame(comp_dat_unval, Y = 1)
      colnames(comp_dat_y0)[length(colnames(comp_dat_y0))] <- colnames(comp_dat_y1)[length(colnames(comp_dat_y1))] <- Y_val
      comp_dat_unval <- data.matrix(cbind(rbind(comp_dat_y0, comp_dat_y1), 
                                          k = rep(rep(seq(1,m), each = (n-n_v)), times = 2)))
      
      comp_dat_all <- rbind(comp_dat_val, comp_dat_unval)
      
      sigma <- sd(VCCC_dat$SQRT_CD4_VAL10)
      theta0 <- matrix(0, nrow = 4, ncol = 1)
      gamma0 <- matrix(0, nrow = (2 + 3*2), ncol = 1)
      lambda0 <- matrix(0, nrow = (6 + 5 + 4), ncol = 1)
      
      mle_res <- nlm(f = observed_data_loglik, p = c(theta0, gamma0, lambda0, sigma), hessian = TRUE, 
                     n = nrow(data), n_v = sum(V), m = m, Y_unval="ADE_WITHIN1YR_UNVAL", Y_val="ADE_WITHIN1YR_VAL", 
                     contX_unval="SQRT_CD4_UNVAL10", contX_val="SQRT_CD4_VAL10", 
                     binX_unval=c("PRIOR_ADE_UNVAL", "ART_NAIVE_ENROLL_UNVAL"), binX_val=c("PRIOR_ADE_VAL", "ART_NAIVE_ENROLL_VAL"), C=NULL, 
                     comp_dat_all = comp_dat_all, print.level = 0, iterlim = 1000
      )
      
      theta_hat <- mle_res$estimate[1:4]
      cov_theta_hat <- tryCatch(expr = sqrt(diag(solve(mle_res$hessian)))[1:4],  error = function(err) {NA})
      
      sim_res$alpha[r] <- theta_hat[1]; sim_res$alphaSE[r] <- cov_theta_hat[1]
      sim_res$betaCD4[r] <- theta_hat[2]; sim_res$betaCD4SE[r] <- cov_theta_hat[2]
      sim_res$betaPriorADE[r] <- theta_hat[3]; sim_res$betaPriorADE[r] <- cov_theta_hat[3]
      sim_res$betaARTNaive[r] <- theta_hat[4]; sim_res$betaARTNaive[r] <- cov_theta_hat[4]
      sim_res$converged[r] <- mle_res$code == 1
      
      write.csv(sim_res, paste0("~sarahlotspeich/Downloads/VCCC_MLE/VCCC_MLE_", nv, "_CaseControl_seed", sim_seed, ".csv"), row.names = FALSE) #
      print(paste("Replication", r, "complete."))
    }
  }
}