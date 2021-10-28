###############################################################
# Simulation setup for Table S7 -------------------------------
# Errors in binary outcome, -----------------------------------
# Continuous covariate (multiplicative) -----------------------
# Varied differential error max -------------------------------
###############################################################

set.seed(918)

# Set sample sizes --------------------------------------------
N <- 1000 # Phase-I = N
n <- 250 # Phase-II/audit size = n

# Generate true values Y, Xb, Xa --------------------------------
Xa <- rbinom(n = N, size = 1, prob = 0.25)
Xb <- rnorm(n = N, mean = 0, sd = 1)
Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(- 1 + Xb - 0.5 * Xa))) ^ (- 1))

# Generate error-prone Xb* = Xb x U -----------------------------
## For U ~ Unif(0, eta0) if C = 0, ~ Unif(0, eta1) if C = 1 ---
eta0 <- 1; eta1 <- 2
U <- runif(n = N, min = 0, max = ifelse(Xa == 0, eta0, eta1))
Xbstar <- Xb * U

# Generate error-prone Y* from error model P(Y*|Xb*,Y,Xb,Xa) ---
sensY <- 0.95; specY <- 0.90
theta0 <- - log(specY / (1 - specY))
theta1 <- - theta0 - log((1 - sensY) / sensY)
Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(- (theta0 + Xbstar + theta1 * Y + Xb - 0.5 * Xa))) ^ (- 1))

# Draw naive case-control audit of size n ---------------------
## V is a TRUE/FALSE vector where TRUE = validated ------------
V <- seq(1, N) %in% c(sample(x = which(Ystar == 0), size = 0.5 * n, replace = FALSE),
                      sample(x = which(Ystar == 1), size = 0.5 * n, replace = FALSE))

# Build dataset ------------------------------------------------
sdat <- cbind(Y, Xb, Ystar, Xbstar, Xa, V)
# Make Phase-II variables Y, Xb NA for unaudited subjects -------
sdat[!V, c("Y", "Xb")] <- NA

# Fit models ---------------------------------------------------
## (1) Naive model ---------------------------------------------
naive <- glm(Ystar ~ Xbstar + Xa, family = "binomial", data = sdat)
beta_naive <- naive$coefficients[2]
se_naive <- sqrt(diag(cov(naive)))[2]

## (2) Complete data -------------------------------------------
cd <- glm(Y[V] ~ Xb[V] + Xa[V], family = "binomial", data = sdat)
beta_cd <- cd$coefficients[2]
se_cd <- sqrt(diag(cov(cd)))[2]

## (3) Horvitz Thompson ------------------------------------
sample_wts <- ifelse(Ystar[V] == 0, 1 / ((0.5 * n) / (table(Ystar)[1])), 1 / ((0.5 * n) / (table(Ystar)[2])))
ht <- glm(Y[V] ~ Xb[V] + Xa[V], family = "binomial", weights = sample_wts)
beta_ht <- ht$coefficients[2]
se_ht <- sqrt(diag(sandwich::sandwich(ht)))[2]

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
sdat <- cbind(id = 1:N, sdat, naive_infl)
library(survey)
sstudy <- twophase(id = list(~id, ~id),
                   data = data.frame(sdat),
                   strat = list(NULL, ~Ystar),
                   subset = ~V)

# Calibrate raking weights to the sum of the naive influence functions ----------------
scal <- calibrate(sstudy, ~ if1 + if2 + if3, phase = 2, calfun = "raking")
# Fit analysis model using calibrated weights -----------------------------------------
rake <- svyglm(Y ~ Xb + Xa, family = "binomial", design = scal)
beta_rake <- rake$coefficients[2]
se_rake <- sqrt(diag(vcov(rake)))[2]

## (5) SMLE ------------------------------------------------
### Construct B-spline basis -------------------------------
### We chose cubic B-splines, with 20 df for N = 1000 ------
nsieve <- 20
B <- matrix(0, nrow = N, ncol = nsieve)
B[which(Xa == 0),1:(0.75 * nsieve)] <- splines::bs(x = Xbstar[which(Xa == 0)], df = 0.75 * nsieve, Boundary.knots = range(Xbstar[which(Xa == 0)]), intercept = TRUE)
B[which(Xa == 1),(0.75 * nsieve + 1):nsieve] <- splines::bs(x = Xbstar[which(Xa == 1)], df = 0.25 * nsieve, Boundary.knots = range(Xbstar[which(Xa == 1)]), intercept = TRUE)
colnames(B) <- paste0("bs", seq(1, nsieve))
sdat <- cbind(sdat, B)

### R package: implementation of proposed SMLE approach ----
### To download the package, run: devtools::install_github("sarahlotspeich/logreg2ph")
library("logreg2ph")
smle <- logreg2ph(Y_unval = "Ystar",
                  Y_val = "Y",
                  X_unval = "Xbstar",
                  X_val = "Xb",
                  C = "Xa",
                  Validated = "V",
                  Bspline = colnames(B),
                  data = sdat,
                  noSE = FALSE,
                  MAX_ITER = 1000,
                  TOL = 1E-4)
beta_smle <- smle$Coefficients$Coefficient[2]
se_smle <- smle$Coefficients$SE[2]
