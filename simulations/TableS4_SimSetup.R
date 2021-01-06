###########################################################
# Simulation setup for Table S4 ---------------------------
# Classical measurement error: Errors in continuous -------
# Covariate with binary error-free covariate --------------
###########################################################

set.seed(918)

# Set sample sizes ----------------------------------------
N <- 1000 # Phase-I = N
n <- 250 # Phase-II/audit size = n

# True parameter values for P(Y|X,Z) ----------------------
beta0 <- - 1
beta1 <- - 2 # vary between [-2, 2]
beta2 <- - 0.50

# Generate true values Y, X, Z ----------------------------
Z <- rnorm(n = N, mean = 0, sd = 1)
X <- rbinom(n = N, size = 1,
            prob = (1 + exp(- (0))) ^ (- 1))
Y <- rbinom(n = N, size = 1,
            prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z))) ^ (- 1))



# Choose audit design: SRS or -----------------------------
## Unvalidated case-control: case-control based on Y^* ----
audit <- "SRS" #or "Validated case-control"

# Draw audit of size n based on design --------------------
## V is a TRUE/FALSE vector where TRUE = validated --------
if(audit == "SRS") {
  V <- seq(1, N) %in% sample(x = seq(1, N), size = n, replace = FALSE)
}
if(audit == "Validated case-control") {
  V <- seq(1, N) %in% c(sample(x = which(Y == 0), size = 0.5 * n, replace = FALSE),
                        sample(x = which(Y == 1), size = 0.5 * n, replace = FALSE))
}

# Build dataset --------------------------------------------
sdat <- cbind(Y, X, Xstar, Z, V)
# Make Phase-II variable X NA for unaudited subjects -------
sdat[!V, c("X")] <- NA

# Fit models -----------------------------------------------
## (1) Naive model -----------------------------------------
naive <- glm(Y ~ Xstar + Z, family = "binomial", data = data.frame(sdat))
beta_naive <- naive$coefficients[2]
se_naive <- sqrt(diag(vcov(naive)))[2]

## (2) Complete case model ---------------------------------
cc <- glm(Y[V] ~ X[V] + Z[V], family = "binomial")
beta_cc <- cc$coefficients[2]
se_cc <- sqrt(diag(vcov(cc)))[2]

## (3) Horvitz-Thompson (HT) estimator ---------------------
## Note: if audit = "SRS", then CC = HT --------------------
if (audit == "Validated case-control") {
  sample_wts <- ifelse(Y[V] == 0, 1 / ((0.5 * n) / (table(Y)[1])), 1 / ((0.5 * n) / (table(Y)[2])))
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
sdat <- cbind(id = 1:N, sdat, naive_infl)
library(survey)
if (audit == "SRS") {
  sstudy <- twophase(id = list(~id, ~id),
                     data = data.frame(sdat),
                     subset = ~V)
} else if (audit == "Validated case-control") {
  sstudy <- twophase(id = list(~id, ~id),
                     data = data.frame(sdat),
                     strat = list(NULL, ~Y),
                     subset = ~V)
}

# Calibrate raking weights to the sum of the naive influence functions ----------------
scal <- calibrate(sstudy, ~ if1 + if2 + if3, phase = 2, calfun = "raking")
# Fit analysis model using calibrated weights -----------------------------------------
rake <- svyglm(Y ~ X + Z, family = "binomial", design = scal)
beta_rake <- rake$coefficients[2]
se_rake <- sqrt(diag(vcov(rake)))[2]

## (5) SMLE ------------------------------------------------
### Construct B-spline basis -------------------------------
### For N = 1000, we select 32 sieves/total df -------------
### Let q = 3, cubic B-spline basis ------------------------
nsieve <- 32
B <- matrix(0, nrow = N, ncol = nsieve)
### Stratify our B-splines on binary C ---------------------
B[which(Z == 0),1:(0.75 * nsieve)] <- splines::bs(x = Xstar[which(Z == 0)], df = 0.75 * nsieve, Boundary.knots = range(Xstar[which(Z == 0)]), intercept = TRUE)
B[which(Z == 1),(0.75 * nsieve + 1):nsieve] <- splines::bs(x = Xstar[which(Z == 1)], df = 0.25 * nsieve, Boundary.knots = range(Xstar[which(Z == 1)]), intercept = TRUE)
colnames(B) <- paste0("bs", seq(1, nsieve))
sdat <- cbind(sdat, B)

### R package: implementation of proposed SMLE approach ----
### To download the package, run: devtools::install_github("sarahlotspeich/logreg2ph")
library("logreg2ph")
smle <- logreg2ph(Y_unval = NULL,
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
