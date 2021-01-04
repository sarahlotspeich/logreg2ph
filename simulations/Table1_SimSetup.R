###########################################################
# Simulation setup for Table 1 ----------------------------
# Errors in binary outcome, binary covariate --------------
###########################################################

set.seed(918)

# Set sample sizes ----------------------------------------
N <- 1000 # Phase-I = N
n <- 250 # Phase-II/audit size = n

# True parameter values for P(Y|X,Z) ----------------------
beta0 <- - 0.65; beta1 <- - 0.20; beta2 <- - 0.10

# Generate true values Y, X, Z ----------------------------
Z <- rbinom(n = N, size = 1,
            prob = 0.25)
X <- rbinom(n = N, size = 1,
            prob = (1 + exp(- (0))) ^ (- 1))
Y <- rbinom(n = N, size = 1,
            prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z))) ^ (- 1))

# Parameters for error model P(X*|X,Z) --------------------
## Set sensitivity/specificity of X* = 0.75 ---------------
sensX <- specX <- 0.75
delta0 <- - log(specX / (1 - specX))
delta1 <- - delta0 - log((1 - sensX) / sensX)
delta2 <- 0 # We assume conditional independence of Y and X* given X
delta3 <- 0.5

# Generate error-prone X* from error model P(X*|X,Z) ------
Xstar <- rbinom(n = N, size = 1,
                 prob = (1 + exp(- (delta0 + delta1 * X + delta2 * Y + delta3 * Z))) ^ (- 1))

# Parameters for error model P(Y*|X*,Y,X,Z) ---------------
## Set sensitivity/specificity of Y* = 0.95, 0.90 ---------
sensY <- 0.95
specY <- 0.90
theta0 <- - log(specY / (1 - specY))
theta1 <- - theta0 - log((1 - sensY) / sensY)
theta2 <- - 0.2
theta3 <- - 0.2
theta4 <- - 0.1

# Generate error-prone Y* from error model P(Y*|X*,Y,X,Z) --
Ystar <- rbinom(n = N, size = 1,
                prob = (1 + exp(- (theta0 + theta1 * Y + theta2 * X + theta3 * Xstar + theta4 * Z))) ^ (- 1))

# Choose audit design: SRS or -----------------------------
## Unvalidated case-control: case-control based on Y^* ----
audit <- "SRS" #or "Unvalidated case-control"

# Draw audit of size n based on design --------------------
## V is a TRUE/FALSE vector where TRUE = validated --------
if(audit == "SRS")
{
    V <- seq(1, N) %in% sample(x = seq(1, N), size = n, replace = FALSE)
}
if(audit == "Unvalidated case-control")
{
    V <- seq(1, N) %in% c(sample(x = which(Ystar == 0), size = 0.5 * n, replace = FALSE),
                          sample(x = which(Ystar == 1), size = 0.5 * n, replace = FALSE))
}

# Build dataset --------------------------------------------
sdat <- cbind(Y, X, Ystar, Xstar, Z, V)
# Make Phase-II variables Y, X NA for unaudited subjects ---
sdat[!V, c("Y", "X")] <- NA

# Fit models -----------------------------------------------
## (1) Naive model -----------------------------------------
naive <- glm(Ystar ~ Xstar + Z, family = "binomial", data = data.frame(sdat))
beta_naive <- naive$coefficients[2]
se_naive <- sqrt(diag(vcov(naive)))[2]

## (2) Complete case model ---------------------------------
cc <- glm(Y[V] ~ X[V] + Z[V], family = "binomial")
beta_cc <- cc$coefficients[2]
se_cc <- sqrt(diag(vcov(cc)))[2]

## (3) Horvitz-Thompson (HT) estimator ---------------------
## Note: if audit = "SRS", then CC = HT --------------------
if (audit == "Unvalidated case-control") {
  library(sandwich)
  sample_wts <- ifelse(Ystar[V] == 0, 1 / ((0.5 * n) / (table(Ystar)[1])), 1 / ((0.5 * n) / (table(Ystar)[2])))
  ht <- glm(Y[V] ~ X[V] + Z[V], family = "binomial",
            weights = sample_wts)
  beta_ht <- ht$coefficients[2]
  se_ht <- sqrt(diag(sandwich(ht)))[2]
}

## (4) Generalized raking ----------------------------------
### Influence function for logistic regression
### Taken from: https://github.com/T0ngChen/multiwave/blob/master/sim.r
inf.fun <- function(fit) {
  dm <- model.matrix(fit)
  Ihat <- (t(dm) %*% (dm * fit$fitted.values * (1 - fit$fitted.values))) / nrow(dm)
  ## influence function
  infl <- (dm * resid(fit, type = "response")) %*% solve(Ihat)
  infl
}

naive_infl <- inf.fun(naive) # error-prone influence functions based on naive model
colnames(naive_infl) <- paste0("if", 1:3)

# Add naive influence functions to sdat -----------------------------------------------
sdat <- cbind(sdat, naive_infl)
library(survey)
if (audit == "SRS") {
  sstudy <- twophase(id = list(~1, ~1),
                     data = data.frame(sdat),
                     subset = ~V)
} else if (audit == "Unvalidated case-control") {
  sstudy <- twophase(id = list(~id, ~id),
                     data = data.frame(sdat),
                     strat = list(NULL, ~Ystar),
                     subset = ~V)
}

# Calibrate raking weights to the sum of the naive influence functions ----------------
scal <- calibrate(sstudy, ~ if1 + if2 + if3, phase = 2, calfun = "raking")
# Fit analysis model using calibrated weights -----------------------------------------
rake <- svyglm(Y ~ X + Z, family = "binomial", design = scal)
beta_rake <- rake$coefficients[2]
se_rake <- sqrt(diag(vcov(rake)))[2]

## (5) MLE -------------------------------------------------
### Script: two-phase log-likelihood specification adapted from Tang et al. (2015)
### Get the script here https://github.com/sarahlotspeich/logreg2ph/blob/master/simulations/Tang_twophase_loglik_binaryX.R
source("Tang_twophase_loglik_binaryX.R")
fit_Tang <- nlm(f = Tang_twophase_loglik,
                p = rep(0, 12),
                hessian = TRUE,
                Val = "V",
                Y_unval = "Ystar",
                Y_val="Y",
                X_unval = "Xstar",
                X_val = "X",
                C = "Z",
                data = sdat)
beta_mle <- fit_Tang$estimate[10]
se_mle <- sqrt(diag(solve(fit_Tang$hessian)))[10]

## (6) SMLE ------------------------------------------------
### Construct B-spline basis -------------------------------
### Since X* and Z are both binary, reduces to indicators --
nsieve <- 4
B <- matrix(0, nrow = N, ncol = nsieve)
B[which(Z == 0 & Xstar == 0), 1] <- 1
B[which(Z == 0 & Xstar == 1), 2] <- 1
B[which(Z == 1 & Xstar == 0), 3] <- 1
B[which(Z == 1 & Xstar == 1), 4] <- 1
colnames(B) <- paste0("bs", seq(1, nsieve))
sdat <- cbind(sdat, B)

### R package: implementation of proposed SMLE approach ----
### To download the package, run: devtools::install_github("sarahlotspeich/logreg2ph")
library("logreg2ph")
smle <- logreg2ph(Y_unval = "Ystar",
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
beta_smle <- smle$Coefficients$Coefficient[2]
se_smle <- smle$Coefficients$SE[2]
