### NOTE: The package is work in progress, and features might be deleted or changed altogether

# Continuous Time Stochastic Modelling for R using Template Model Builder (ctsmrTMB)

`ctsmrTMB` is an R package for parameter estimation, state filtration and forecasting in stochastic state space models, heavily inspired by [Continuous Time Stochastic Modelling](https://ctsm.info). 
The package is a user-friendly wrapper for [Template Model Builder](https://github.com/kaskr/adcomp) that frees the user from writing
the required C++ file containing the (negative log) likelihood function themselves. Instead, the C++ script is generated automatically based on a model specified by the user using the provided R6 `ctsmrTMB` class object.


