###########################################################
# Simulation setup for Table S2 ---------------------------
# Errors in binary outcome, -------------------------------
# Continuous covariate (additive, nonzero mean) -----------
# Varied differential error means -------------------------
###########################################################

set.seed(918)

# Set sample sizes ----------------------------------------
N <- 1000 # Phase-I = N
n <- 250 # Phase-II/audit size = n

# True parameter values for P(Y|X,Z) ----------------------
beta0 <- -1; beta1 <- 1; beta2 <- -0.5

# Generate true values Y, X, Z ----------------------------
Z <- rbinom(n = N, size = 1, 
            prob = 0.25)
X <- rbinom(n = N, size = 1, 
            prob = (1 + exp(- (0))) ^ (- 1))
Y <- rbinom(n = N, size = 1, 
            prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z))) ^ (- 1))

# Generate error-prone X* = X + U -------------------------
## For U ~ N(mu0, var) if C = 0, ~ N(mu1, var) if C = 1 ---
mu0 <- 0; mu1 <- 1 
s2U <- 0.1 
U <- rnorm(n = N, mean = ifelse(Z == 0, mu0, mu1), sd = sqrt(sU))
Xstar <- X + U

# Parameters for error model P(Y*|X*,Y,X,Z) ---------------
## Set sensitivity/specificity of Y* ----------------------
sensY <- 0.95; specY <- 0.90
theta0 <- - log(specY / (1 - specY))
theta1 <- - theta0 - log((1 - sensY) / sensY)
theta2 <- 1; theta3 <- 1; theta4 <- - 0.5

# Generate error-prone Y* from error model P(Y*|X*,Y,X,Z) --
Ystar <- rbinom(n = N, size = 1, 
                prob = (1 + exp(- (theta0 + theta1 * Y + theta2 * X + theta3 * Xstar + theta4 * Z))) ^ (- 1))

# Choose audit design: ------------------------------------
## Unvalidated case-control: case-control based on Y^* ----
audit <- "Unvalidated case-control"

# Draw audit of size n based on design --------------------
## V is a TRUE/FALSE vector where TRUE = validated --------
V <- seq(1, N) %in% c(sample(x = which(Ystar == 0), size = 0.5 * n, replace = FALSE), 
                      sample(x = which(Ystar == 1), size = 0.5 * n, replace = FALSE))

# Build dataset --------------------------------------------
sdat <- cbind(Y, X, Ystar, Xstar, Z, V)
# Make Phase-II variables Y, X NA for unaudited subjects ---
sdat[!V, c("Y", "X")] <- NA

# Fit models -----------------------------------------------
## (1) Naive model -----------------------------------------
naive <- glm(Ystar ~ Xstar + Z, family = "binomial", data = sdat)

## (2) Horvitz Thompson ------------------------------------
ht <- glm(Y[V] ~ X[V] + Z[V], family = "binomial", data = sdat,
          weights = ifelse(Ystar[V] == 0, 
                           1 / ((0.5 * n) / (table(Ystar)[1])), 
                           1 / ((0.5 * n) / (table(Ystar)[2]))))
beta_ht <- ht$coefficients[2]
se_ht <- sqrt(diag(sandwich::sandwich(ht)))[2] 

## (3) SMLE ------------------------------------------------
### Construct B-spline basis -------------------------------
### We chose cubic B-splines, with 20 df for N = 1000 ------
nsieve <- 20
B <- matrix(0, nrow = N, ncol = nsieve)
B[which(Z == 0),1:(0.75 * nsieve)] <- splines::bs(x = Xstar[which(Z == 0)], df = 0.75 * nsieve, Boundary.knots = range(Xstar[which(Z == 0)]), intercept = TRUE)
B[which(Z == 1),(0.75 * nsieve + 1):nsieve] <- splines::bs(x = Xstar[which(Z == 1)], df = 0.25 * nsieve, Boundary.knots = range(Xstar[which(Z == 1)]), intercept = TRUE)
colnames(B) <- paste0("bs", seq(1, nsieve))
sdat <- cbind(sdat, B)

### Script: implementation of proposed SMLE approach -------
### Get the script here https://github.com/sarahlotspeich/logreg_2ph/blob/master/logreg_2ph.R
source("logreg_2ph.R")

smle <- logreg_2ph(Y_unval = "Ystar", Y_val = "Y", X_unval = "Xstar", X_val = "X", C = "Z", Validated = "V", Bspline = colnames(B), 
                   data = sdat, noSE = FALSE, MAX_ITER = 1000, rescale = FALSE, TOL = 1E-4)
beta_smle <- smle$Coefficients$Coefficient[2]
se_smle <- smle$Coefficients$SE[2]