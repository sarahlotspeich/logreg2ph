context("logreg2ph complete")

library(testthat)
library(logreg2ph)
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)


## INITIALIZE FUNCTIONS
if (!exists("logreg2ph"))
{
	source("R/hessian_row.R")
	source("R/logreg2ph.R")
	source("R/observed_data_loglik.R")
	source("R/pl_theta.R")
	source("R/profile_out.R")
	Rcpp::sourceCpp("src/fast_functions.cpp")
}



## INITIALIZE ALL VARIABLES FOR TABLE 1 TESTS
set.seed(918)

# Set sample sizes ----------------------------------------
N <- 1000 # Phase-I = N
n <- 250 # Phase-II/audit size = n

# Generate true values Y, Xb, Xa --------------------------
Xa <- rbinom(n = N, size = 1, prob = 0.25)
Xb <- rbinom(n = N, size = 1, prob = 0.5)
Y <- rbinom(n = N, size = 1,prob = (1 + exp(-(- 0.65 - 0.2 * Xb - 0.1 * Xa))) ^ (- 1))

# Generate error-prone Xb* from error model P(Xb*|Xb,Xa) --
sensX <- specX <- 0.75
delta0 <- - log(specX / (1 - specX))
delta1 <- - delta0 - log((1 - sensX) / sensX)
Xbstar <- rbinom(n = N, size = 1,
	prob = (1 + exp(- (delta0 + delta1 * Xb + 0.5 * Xa))) ^ (- 1))

# Generate error-prone Y* from error model P(Y*|Xb*,Y,Xb,Xa)
sensY <- 0.95
specY <- 0.90
theta0 <- - log(specY / (1 - specY))
theta1 <- - theta0 - log((1 - sensY) / sensY)
Ystar <- rbinom(n = N, size = 1,
	prob = (1 + exp(- (theta0 - 0.2 * Xbstar + theta1 * Y - 0.2 * Xb - 0.1 * Xa))) ^ (- 1))

# Choose audit design: SRS or -----------------------------
## Naive case-control: case-control based on Y^* ----------
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




test_that("Table 1", {

	# Build dataset --------------------------------------------
	sdat <- cbind(Y, Xb, Ystar, Xbstar, Xa, V)
	# Make Phase-II variables Y, Xb NA for unaudited subjects ---
	sdat[!V, c("Y", "Xb")] <- NA


	## (6) SMLE ------------------------------------------------
	### Construct B-spline basis -------------------------------
	### Since Xb* and Xa are both binary, reduces to indicators --
	nsieve <- 4
	B <- matrix(0, nrow = N, ncol = nsieve)
	B[which(Xa == 0 & Xbstar == 0), 1] <- 1
	B[which(Xa == 0 & Xbstar == 1), 2] <- 1
	B[which(Xa == 1 & Xbstar == 0), 3] <- 1
	B[which(Xa == 1 & Xbstar == 1), 4] <- 1
	colnames(B) <- paste0("bs", seq(1, nsieve))
	sdat <- cbind(sdat, B)

	expect_message(smle <- logreg2ph(Y_unval = "Ystar",
		Y_val = "Y",
		X_unval = "Xbstar",
		X_val = "Xb",
		C = "Xa",
		Validated = "V",
		Bspline = colnames(B),
		data = sdat,
		noSE = FALSE,
		MAX_ITER = 1000,
		TOL = 1E-4), "Analysis model assumed main effects only")

	# July 16 2021: 0.78 sec

	# expect_lt(max(abs(smle$Coefficients - matrix(c(-0.587801556603326, -0.231437717735961, 0.139139424615862, 0.159867158182001, 0.240522213901066, 0.203146591074279), nrow=3))), 1e-9)
	# expect_lt(max(abs(smle$I_theta - c(-125.250617048664, -63.4359041852122, -32.8821658440575, -63.4359041852122, -49.4155484420844, -16.8277466843847, -32.8821658440575, -16.8277466843847, -32.8658701123459))), 1e-9)
	# expect_true(smle$converged)
	# expect_true(smle$se_converged)
	# expect_equal(smle$converged_msg, "Converged")
	# expect_equal(smle$iterations, 42)
	# expect_equal(smle$od_loglik_at_conv, -847.635948985586)

	expect_equal(smle, logreg2ph::logreg2ph(Y_unval = "Ystar",
	                             Y_val = "Y",
	                             X_unval = "Xbstar",
	                             X_val = "Xb",
	                             C = "Xa",
	                             Validated = "V",
	                             Bspline = colnames(B),
	                             data = sdat,
	                             noSE = FALSE,
	                             MAX_ITER = 1000,
	                             TOL = 1E-4))
})


test_that("Table 1 - no errors Y", {

	# Build dataset --------------------------------------------
	sdat <- cbind(Y, Xb, Xbstar, Xa, V)
	# Make Phase-II variables Xb NA for unaudited subjects ---
	sdat[!V, c("Xb")] <- NA


	## (6) SMLE ------------------------------------------------
	### Construct B-spline basis -------------------------------
	### Since Xb* and Xa are both binary, reduces to indicators --
	nsieve <- 4
	B <- matrix(0, nrow = N, ncol = nsieve)
	### Stratify our B-splines on binary C ---------------------
	B[which(Xa == 0 & Xbstar == 0), 1] <- 1
	B[which(Xa == 0 & Xbstar == 1), 2] <- 1
	B[which(Xa == 1 & Xbstar == 0), 3] <- 1
	B[which(Xa == 1 & Xbstar == 1), 4] <- 1
	colnames(B) <- paste0("bs", seq(1, nsieve))
	sdat <- cbind(sdat, B)

	smle <- logreg2ph(Y_unval = NULL,
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

	expect_equal(smle, logreg2ph::logreg2ph(Y_unval = NULL,
	                             Y_val = "Y",
	                             X_unval = "Xbstar",
	                             X_val = "Xb",
	                             C = "Xa",
	                             Validated = "V",
	                             Bspline = colnames(B),
	                             data = sdat,
	                             noSE = FALSE,
	                             MAX_ITER = 1000,
	                             TOL = 1E-4))
})



test_that("Table 1 - no errors X", {

	# Build dataset --------------------------------------------
	sdat <- cbind(Y, Xb, Ystar, Xa, V)
	# Make Phase-II variables Y NA for unaudited subjects ---
	sdat[!V, c("Y")] <- NA


	## (6) SMLE ------------------------------------------------
	### Construct B-spline basis -------------------------------
	### Since Xb* and Xa are both binary, reduces to indicators --
	nsieve <- 4
	B <- matrix(0, nrow = N, ncol = nsieve)
	### Stratify our B-splines on binary C ---------------------
	B[which(Xa == 0 & Ystar == 0), 1] <- 1
	B[which(Xa == 0 & Ystar == 1), 2] <- 1
	B[which(Xa == 1 & Ystar == 0), 3] <- 1
	B[which(Xa == 1 & Ystar == 1), 4] <- 1
	colnames(B) <- paste0("bs", seq(1, nsieve))
	sdat <- cbind(sdat, B)

	smle <- logreg2ph(Y_unval = "Ystar",
	                  Y_val = "Y",
	                  X_unval = NULL,
	                  X_val = "Xb",
	                  C = "Xa",
	                  Validated = "V",
	                  Bspline = colnames(B),
	                  data = sdat,
	                  noSE = FALSE,
	                  MAX_ITER = 1000,
	                  TOL = 1E-4)

	expect_equal(smle, logreg2ph::logreg2ph(Y_unval = "Ystar",
	                  Y_val = "Y",
	                  X_unval = NULL,
	                  X_val = "Xb",
	                  C = "Xa",
	                  Validated = "V",
	                  Bspline = colnames(B),
	                  data = sdat,
	                  noSE = FALSE,
	                  MAX_ITER = 1000,
	                  TOL = 1E-4))
})



test_that("Table 4", {
	skip_on_cran()
  	# this test can take several minutes to run
  	skip("skipping slow test, comment out this line if you want to run it")

	###########################################################
	# Simulation setup for Table 4 ----------------------------
	# Errors in binary outcome/continuous covariate (additive)-
	# Varied outcome error rates ------------------------------
	###########################################################

	set.seed(918)

	# Set sample sizes ----------------------------------------
	N <- 1000 # Phase-I = N
	n <- 250 # Phase-II/audit size = n

	# Generate true values Y, Xb, Xa ----------------------------
	Xa <- rbinom(n = N, size = 1, prob = 0.25)
	Xb <- rnorm(n = N, mean = 0, sd = 1)
	Y <- rbinom(n = N, size = 1, prob = (1 + exp(-(- 1 + Xb - 0.5 * Xa))) ^ (- 1))

	# Generate error-prone Xb* = Xb + U -------------------------
	## For U ~ N(mean, var) -----------------------------------
	muU <- 0 ## mean muU = 0 for unbiased, != 0 for biased ----
	s2U <- 0.1 ## variance ------------------------------------
	U <- rnorm(n = N, mean = muU, sd = sqrt(s2U))
	Xbstar <- Xb + U

	# Parameters for error model P(Y*|Xb*,Y,Xb,Xa) ---------------
	## Varied outcome error rates,
	## assumed specificity = sensitivity - 0.05
	sensY <- 0.95 # 0.95, 0.85, 0.75, 0.65, 0.55
	specY <- 0.90 # 0.90, 0.80, 0.70, 0.60, 0.50
	theta0 <- - log(specY / (1 - specY))
	theta1 <- - theta0 - log((1 - sensY) / sensY)
	theta2 <- 1; theta3 <- 1; theta4 <- - 0.5

	# Generate error-prone Y* from error model P(Y*|Xb*,Y,Xb,Xa) --
	Ystar <- rbinom(n = N, size = 1,
	                prob = (1 + exp(- (theta0 + theta1 * Y + theta2 * Xb + theta3 * Xbstar + theta4 * Xa))) ^ (- 1))

	# Choose audit design: SRS or -----------------------------
	## Unvalidated case-control: case-control based on Y^* ----
	## Table 2 was SRS, Table 3 was Unvalidated case-control --
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
	sdat <- data.frame(cbind(Y, Xb, Ystar, Xbstar, Xa, V))
	# Make Phase-II variables Y, Xb NA for unaudited subjects ---
	sdat[!V, c("Y", "Xb")] <- NA

	## (5) SMLE ------------------------------------------------
	### Construct B-spline basis -------------------------------
	### We chose cubic B-splines, ------------------------------
	### with 20 df for N = 1000 and 24 df for N = 2000 ---------
	nsieve <- 20
	B <- matrix(0, nrow = N, ncol = nsieve)
	B[which(Xa == 0),1:(0.75 * nsieve)] <- splines::bs(x = Xbstar[which(Xa == 0)], df = 0.75 * nsieve, Boundary.knots = range(Xbstar[which(Xa == 0)]), intercept = TRUE)
	B[which(Xa == 1),(0.75 * nsieve + 1):nsieve] <- splines::bs(x = Xbstar[which(Xa == 1)], df = 0.25 * nsieve, Boundary.knots = range(Xbstar[which(Xa == 1)]), intercept = TRUE)
	colnames(B) <- paste0("bs", seq(1, nsieve))
	sdat <- cbind(sdat, B)

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

	expect_equal(smle, logreg2ph::logreg2ph(Y_unval = "Ystar",
	                  Y_val = "Y",
	                  X_unval = "Xbstar",
	                  X_val = "Xb",
	                  C = "Xa",
	                  Validated = "V",
	                  Bspline = colnames(B),
	                  data = sdat,
	                  noSE = FALSE,
	                  MAX_ITER = 1000,
	                  TOL = 1E-4))


	# without cpp:  311.38 sec
	# June 24 2021: 220.46 sec
	# June 29 2021: 176.23 sec
	# July 15 2021: 155.47 sec
	# July 22 2021: 141.77 sec
	# July 23 2021: 114.31 sec
	# July 26 2021: 118.08 sec
	#       min       lq     mean   median       uq      max  neval
	#   97.01602 102.7917 107.2892 106.3271 111.1261 118.6213    10
})


test_that("Table S4 failing", {
	# Build dataset --------------------------------------------
	sdat <- cbind(Y, Xb, Xbstar, Xa, V)
	# Make Phase-II variable Xb NA for unaudited subjects -------
	sdat[!V, c("Xb")] <- NA

	## (5) SMLE ------------------------------------------------
	### Construct B-spline basis -------------------------------
	### For N = 1000, we select 32 sieves/total df -------------
	### Let q = 3, cubic B-spline basis ------------------------
	nsieve <- 32
	B <- matrix(0, nrow = N, ncol = nsieve)
	### Stratify our B-splines on binary C ---------------------
	B[which(Xa == 0),1:(0.75 * nsieve)] <- splines::bs(x = Xbstar[which(Xa == 0)], df = 0.75 * nsieve, Boundary.knots = range(Xbstar[which(Xa == 0)]), intercept = TRUE)
	B[which(Xa == 1),(0.75 * nsieve + 1):nsieve] <- splines::bs(x = Xbstar[which(Xa == 1)], df = 0.25 * nsieve, Boundary.knots = range(Xbstar[which(Xa == 1)]), intercept = TRUE)
	colnames(B) <- paste0("bs", seq(1, nsieve))
	sdat <- cbind(sdat, B)

	expect_warning(logreg2ph(Y_unval = NULL,
	                  Y_val = "Y",
	                  X_unval = "Xbstar",
	                  X_val = "Xb",
	                  C = "Xa",
	                  Validated = "V",
	                  Bspline = colnames(B),
	                  data = sdat,
	                  noSE = FALSE,
	                  MAX_ITER = 1000,
	                  TOL = 1E-4), "Empty sieve in validated data. Reconstruct B-spline basis and try again.")

	expect_warning(logreg2ph::logreg2ph(Y_unval = NULL,
	                  Y_val = "Y",
	                  X_unval = "Xbstar",
	                  X_val = "Xb",
	                  C = "Xa",
	                  Validated = "V",
	                  Bspline = colnames(B),
	                  data = sdat,
	                  noSE = FALSE,
	                  MAX_ITER = 1000,
	                  TOL = 1E-4), "Empty sieve in validated data. Reconstruct B-spline basis and try again.")

})
