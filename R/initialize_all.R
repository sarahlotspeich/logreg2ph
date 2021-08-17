# meant for use in development and working with the R notebook

library(tictoc)
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
library(logreg2ph)
library(testthat)

source("R/hessian_row.R")
source("R/logreg2ph.R")
source("R/observed_data_loglik.R")
source("R/pl_theta.R")
source("R/profile_out.R")
Rcpp::sourceCpp("src/fast_functions.cpp")
