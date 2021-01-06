# Efficient odds ratio estimation under two-phase sampling using error-prone data from a multi-national HIV research cohort
## Lotspeich, Shepherd, Amorim, Shaw, and Tao
The complete R package `logreg2ph` and code for the simulation settings included in the paper. 

### Install
To install the package, run the following in your `R` console: `devtools::install_github("sarahlotspeich/logreg2ph")`.

### Simulation settings 
Inside the `simulations` subdirectory, you will find the following: 

  - `Table1_SimSetup.R`: simulations with outcome misclassification and a binary error-prone covariate, intended to inspect increasing Phase I and Phase II sample sizes
  - `Tables2&3_SimSetup.R`: Simulations with outcome misclassification and a continuous covariate with additive errors, intended to inspect varied error variance 
  - `Tables4_SimSetup.R`: Simulations with outcome misclassification and a continuous covariate with additive errors, intended to inspect varied outcome error rates
