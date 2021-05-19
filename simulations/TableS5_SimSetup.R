###########################################################
# Simulation setup for Table S5 ---------------------------
# Behavior of the naive estimator -------------------------
###########################################################

library(magrittr)

set.seed(918)

# Set sample sizes ----------------------------------------
N <- 1000 # Phase-I = N
n <- 250 # Phase-II/audit size = n

# True parameter values for P(Y|X,Z) ----------------------
beta0 <- - 1; beta1 <- 1; beta2 <- - 0.5

# Generate true values Y, X, Z ----------------------------
Z <- rbinom(n = N, size = 1, prob = 0.25)
X <- rnorm(n = N, mean = 0, sd = 1)
Y <- rbinom(n = N, size = 1,
            prob = (1 + exp(- (beta0 + beta1 * X + beta2 * Z)))^(-1))

# Generate error-prone X* = X + U -------------------------
## For U ~ N(muU, s2U) ------------------------------------
muU <- 0
s2U <- 0.1 # 0.1, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2
U <- rnorm(n = N, mean = muU, sd = sqrt(s2U))
Xstar <- X + U

# Setting 1: Errors in X* only ============================
naive_Xonly <- glm(Y ~ Xstar + Z, family = "binomial")

# Setting 2: Errors in X* and Y* ==========================
sensY <- 0.95; specY <- 0.90
gamma0 <- - log(specY / (1 - specY))
gamma1 <- - gamma0 - log((1 - sensY) / sensY)
gamma2 <- gamma3 <- 1
gamma4 <- - 0.5
Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(- (gamma0 + gamma2 * Xstar + gamma1 * Y + gamma3 * X + gamma4 * Z)))^(-1))
naive_Both <- glm(Ystar ~ Xstar + Z, family = "binomial")
