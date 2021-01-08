library(magrittr)

sim_seed <- 918

sensY <- 0.95
specY <- 0.90

prevZ <- 0.25
num_sims <- 1000

beta0 <- - 1
beta1 <- 1
beta2 <- - 0.5

muU <- 0

theta0 <- -log(specY / (1 - specY))
theta1 <- -theta0 - log((1 - sensY) / sensY)
theta2 <- 1
theta3 <- 1
theta4 <- - 0.5

res <- data.frame()

for (sU in c(0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2))
{
    for (N in c(1000))
    {
        try_pv <- c(0.25)
        for (pv in try_pv)
        {
            n <- N*pv
            set.seed(sim_seed)
            
            for(s in 1:num_sims)
            {
                Z <- rbinom(n = N, size = 1, prob = prevZ)
                X <- rnorm(n = N, mean = 0, sd = 1)
                Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(beta0 + beta1*X + beta2*Z)))^(-1))
                U <- rnorm(n = N, mean = muU, sd = sqrt(sU))
                Xstar <- X + U
                
                naive <- glm(Y ~ Xstar + Z, family = "binomial")
                res %<>% dplyr::bind_rows(
                    data.frame(setting = "X only", sigma2U = sU, 
                               beta0 = naive$coefficients[1], beta1 = naive$coefficients[2], beta2 = naive$coefficients[3],
                               fprY = NA, fnrY = NA)
                )
                
                Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(-(theta0 + theta1*Y + theta2*X + theta3*Xstar + theta4*Z)))^(-1))
                naive <- glm(Ystar ~ Xstar + Z, family = "binomial")
                res %<>% dplyr::bind_rows(
                    data.frame(setting = "Y/X both", sigma2U = sU, 
                               beta0 = naive$coefficients[1], beta1 = naive$coefficients[2], beta2 = naive$coefficients[3],
                               fprY = mean(Ystar == 1 & Y == 0), fnrY = mean(Ystar == 0 & Y == 1))
                )
                
                Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(-(theta0 + theta1*Y + theta2*X + theta4*Z)))^(-1))
                naive <- glm(Ystar ~ X + Z, family = "binomial")
                res %<>% dplyr::bind_rows(
                    data.frame(setting = "Y only", sigma2U = sU, 
                               beta0 = naive$coefficients[1], beta1 = naive$coefficients[2], beta2 = naive$coefficients[3],
                               fprY = mean(Ystar == 1 & Y == 0), fnrY = mean(Ystar == 0 & Y == 1))
                )
            }
        }
    }
}