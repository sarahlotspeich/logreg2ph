# Efficient odds ratio estimation under two-phase sampling using error-prone data from a multi-national HIV research cohort
## Lotspeich, Shepherd, Amorim, Shaw, and Tao
The complete R package `logreg2ph` and code for the simulation settings included in the paper. 

### Install
To install the package, run the following in your `R` console: 

```{r}
devtools::install_github("sarahlotspeich/logreg2ph")
```

### Simulation settings 
Inside the `simulations` subdirectory, you will find the following: 

  - `Table1_SimSetup.R`: simulations with outcome misclassification and a binary error-prone covariate, intended to inspect increasing Phase I and Phase II sample sizes
  - `Tables2&3_SimSetup.R`: Simulations with outcome misclassification and a continuous covariate with additive errors, intended to inspect varied error variance 
  - `Table4_SimSetup.R`: Simulations with outcome misclassification and a continuous covariate with additive errors, intended to inspect varied outcome error rates
  - `TableS2_SimSetup.R`: Simulations with outcome misclassification and binary error-prone covariate, intended to inspect robustness of the fully-parametric MLE and proposed SMLE under complex covariate error
  - `TableS3_SimSetup.R`: Simulations with outcome misclassification and binary error-prone covariate, intended to inspect robustness of the fully-parametric MLE and proposed SMLE under complex outcome and covariate error
  - `TableS4_SimSetup.R`: Simulations with continuous covariate with additive errors (classical measurement error), intended to compare performance of proposed SMLE to traditional regression calibration (RC) approach
  - `TableS5_SimSetup.R`: Simulation results for the naive estimator under outcome misclassification and a continuous covariate under varied parameterizations
  - `TableS6_SimSetup.R`: Simulations with outcome misclassification and continuous covariate with differentially biased additive errors, intended to demonstrate robustness of the proposed SMLE to different covariate error types
  - `TableS7_SimSetup.R`: Simulations with outcome misclassification and continuous covariate with differential multiplicative errors, intended to demonstrate robustness of the proposed SMLE to different covariate error types

### Generating Data

Inside each of the files above, you will find code to generate the appropriate data for that simulation setting, e.g., 

```{r, eval = F, tidy = TRUE}
set.seed(918)

# Set sample sizes ----------------------------------------
N <- 1000 # Phase I size = N
n <- 250 # Phase II/audit size = n

# Generate true values Y, Xb, Xa ----------------------------
Xa <- rbinom(n = N, size = 1, prob = 0.25)
Xb <- rbinom(n = N, size = 1, prob = 0.5)
Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(- 0.65 - 0.2 * Xb - 0.1 * Xa))) ^ (- 1))

# Generate error-prone Xb* from error model P(Xb*|Xb,Xa) ------
sensX <- specX <- 0.75
delta0 <- - log(specX / (1 - specX))
delta1 <- - delta0 - log((1 - sensX) / sensX)
Xbstar <- rbinom(n = N, size = 1, prob = (1 + exp(- (delta0 + delta1 * Xb + 0.5 * Xa))) ^ (- 1))

# Generate error-prone Y* from error model P(Y*|Xb*,Y,Xb,Xa) ---
sensY <- 0.95
specY <- 0.90
theta0 <- - log(specY / (1 - specY))
theta1 <- - theta0 - log((1 - sensY) / sensY)
Ystar <- rbinom(n = N, size = 1, prob = (1 + exp(- (theta0 - 0.2 * Xbstar + theta1 * Y - 0.2 * Xb - 0.1 * Xa))) ^ (- 1))
```

Then, the user has the option of two audit designs: simple random sampling (SRS) or 1:1 case-control sampling based on Y* (naive case-control). Based on these designs, the validation indicators V are generated as follows: 

```{r, eval = FALSE}
# Choose audit design: SRS or -----------------------------
## Naive case-control: case-control based on Y^* ----
audit <- "SRS" #or "Naive case-control"

# Draw audit of size n based on design --------------------
## V is a TRUE/FALSE vector where TRUE = validated --------
if(audit == "SRS")
{
    V <- seq(1, N) %in% sample(x = seq(1, N), size = n, replace = FALSE)
}
if(audit == "Naive case-control")
{
    V <- seq(1, N) %in% c(sample(x = which(Ystar == 0), size = 0.5 * n, replace = FALSE),
                          sample(x = which(Ystar == 1), size = 0.5 * n, replace = FALSE))
}
```

Finally, combine the generated data and validation indicators into an analytical dataset: 

```{r, eval = F, tidy = T}
# Build dataset --------------------------------------------
sdat <- cbind(Y, X, Ystar, Xstar, Z, V)
# Make Phase-II variables Y, X NA for unaudited subjects ---
sdat[!V, c("Y", "X")] <- NA
```

### Running Estimator Code

The `R` scripts each contain implementations for the estimators discussed in the paper. Examples of each are provided below:

#### 1. Naive Analysis

```{r, eval = F, tidy = T}
naive <- glm(Ystar ~ Xstar + Z, family = "binomial", data = data.frame(sdat))
beta_naive <- naive$coefficients[2]
se_naive <- sqrt(diag(vcov(naive)))[2]
```

#### 2. Complete-Case Analysis

```{r, eval = F, tidy = T}
cc <- glm(Y[V] ~ X[V] + Z[V], family = "binomial")
beta_cc <- cc$coefficients[2]
se_cc <- sqrt(diag(vcov(cc)))[2]
```

#### 3. Horvitz--Thompson Estimator (for Naive Case-Control Audit Only)

```{r, eval = F, tidy = T}
library(sandwich)
if (audit == "Naive case-control") {
  sample_wts <- ifelse(Ystar[V] == 0, 1 / ((0.5 * n) / (table(Ystar)[1])), 1 / ((0.5 * n) / (table(Ystar)[2])))
  ht <- glm(Y[V] ~ X[V] + Z[V], family = "binomial",
            weights = sample_wts)
  beta_ht <- ht$coefficients[2]
  se_ht <- sqrt(diag(sandwich(ht)))[2]
}
```

#### 4. Generalized Raking Estimator

```{r, eval = F, tidy = T}
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

# Calibrate raking weights to the sum of the naive influence functions ----------------
library(survey)
if (audit == "SRS") {
  sstudy <- twophase(id = list(~id, ~id),
                     data = data.frame(sdat),
                     subset = ~V)
} else if (audit == "Naive case-control") {
  sstudy <- twophase(id = list(~id, ~id),
                     data = data.frame(sdat),
                     strat = list(NULL, ~Ystar),
                     subset = ~V)
}
scal <- calibrate(sstudy, ~ if1 + if2 + if3, phase = 2, calfun = "raking")

# Fit analysis model using calibrated weights -----------------------------------------
rake <- svyglm(Y ~ X + Z, family = "binomial", design = scal)
beta_rake <- rake$coefficients[2]
se_rake <- sqrt(diag(vcov(rake)))[2]
```

#### 5. Maximum Likelihood Estimator (MLE) (for Binary $X^*$ Only)

```{r, eval = F, tidy = T}
# Script: two-phase log-likelihood specification adapted from Tang et al. (2015) named ~/code/Tang_twophase_loglik_binaryX.R
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
```

#### 6. Sieve Maximum Likelihood Estimator (SMLE)

```{r, eval = F, tidy = T}
# Construct B-spline basis -------------------------------
nsieve <- 4
B <- matrix(0, nrow = N, ncol = nsieve)
B[which(Z == 0 & Xstar == 0), 1] <- 1
B[which(Z == 0 & Xstar == 1), 2] <- 1
B[which(Z == 1 & Xstar == 0), 3] <- 1
B[which(Z == 1 & Xstar == 1), 4] <- 1
colnames(B) <- paste0("bs", seq(1, nsieve))
sdat <- cbind(sdat, B)

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
```
