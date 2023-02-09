### NOTE: The package is work in progress, and features might be deleted or changed altogether

# Stochastic Differential Equation Modelling (SDEM) for Template Model Builder (TMB)

`sdemTMB` is an R package for parameter estimation, state filtration and forecasting in stochastic state space models, inspired by [Continuous Time Stochastic Modelling](https://ctsm.info). 
The package is a user-friendly wrapper for [Template Model Builder](https://github.com/kaskr/adcomp) that frees the user from writing
the required C++ file containing the (negative log) likelihood function themselves. Instead, the C++ script is generated automatically based on a model specified by the user using the provided R6 `sdemTMB` class object.  

The package implements the following methods  
1. The `TMB`-style approach where latent states are considered random effects (see e.g. [this example]( https://github.com/kaskr/adcomp/blob/master/tmb_examples/sde_linear.cpp))
2. The (Continous-Discrete) Extended Kalman Filter 
3. The (Continous-Discrete) Unscented Kalman Filter

The main advantage of the kalman filter implementations is a large increase in computational speed (x20), and access to the fixed effects hessian. A distinct advantage of the `TMB`-style implementation is that it allows for non-Gaussian observation noise, but this functionality is not yet implemented in the package.

## Installation

``` r
remotes::install_github(repo="phillipbvetter/sdemTMB")
```
