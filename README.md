# Efficient odds ratio estimation under two-phase sampling using error-prone data from a multi-national HIV research cohort
## Lotspeich, Shepherd, Amorim, Shaw, and Tao
The complete R package `logreg2ph` and code for the simulation settings included in the paper. 

### Install
To install the package, run the following in your `R` console: `devtools::install_github("sarahlotspeich/logreg2ph")`.

### Simulation settings 
Inside the `simulations` subdirectory, you will find the following: 

  - `Table1_SimSetup.R`: simulations with outcome misclassification and a binary error-prone covariate, intended to inspect increasing Phase I and Phase II sample sizes
  - `Tables2&3_SimSetup.R`: Simulations with outcome misclassification and a continuous covariate with additive errors, intended to inspect varied error variance 
  - `Table4_SimSetup.R`: Simulations with outcome misclassification and a continuous covariate with additive errors, intended to inspect varied outcome error rates
  - `TableS2_SimSetup.R`: Simulations with outcome misclassification and binary error-prone covariate, intended to compare performance of the fully-parametric MLE and proposed SMLE under complex covariate error
  - `TableS4_SimSetup.R`: Simulations with continuous covariate with additive errors (classical measurement error), intended to compare performance of proposed SMLE to traditional regression calibration approach
  - `TableS6_SimSetup.R`: Simulations with outcome misclassification and continuous covariate with differentially biased additive errors, intended to demonstrate robustness of the proposed SMLE to different covariate error types
  - `TableS7_SimSetup.R`: Simulations with outcome misclassification and continuous covariate with differential multiplicative errors, intended to demonstrate robustness of the proposed SMLE to different covariate error types
