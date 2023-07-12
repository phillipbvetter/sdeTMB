# Continuous Time Stochastic Modelling for R using Template Model Builder (ctsmrTMB)

`ctsmrTMB` is an R package for parameter estimation, state filtration and forecasting in stochastic state space models, heavily inspired by [Continuous Time Stochastic Modelling](https://ctsm.info). 
The package is a user-friendly wrapper for [Template Model Builder](https://github.com/kaskr/adcomp) that frees the user from writing
the required C++ file containing the (negative log) likelihood function themselves. Instead, the C++ script is generated automatically based on a model specified by the user using the provided R6 `ctsmrTMB` class object.  

The package implements the following methods 
 
1. The `TMB`-style approach where latent states are considered random effects (see e.g. [this example]( https://github.com/kaskr/adcomp/blob/master/tmb_examples/sde_linear.cpp))

2. The (Continous-Discrete) Extended Kalman Filter 

<!--
3. The (Continous-Discrete) Unscented Kalman Filter
-->

The main advantage of the Kalman Filter implementation is a large increase in the computation speed, and access to the fixed effects hessian for improved convergence of the optimization. In these cases TMB just provides automatic differentiation.

A district advantage of the `TMB`-style implementation is its use of the Laplace approximation for likelihood calculations which allows state space formulations where the density of the observation residuals are non-Gaussian. This feature will be implemented in the future.

The package is currently mostly tailored towards the Kalman Filter, with features such as `predict` method for simulations and k-step-ahead predictions without data updates, and `plot` of the fit object from `estimate` to display a basic residuals analysis. These features will eventually be available for the TMB method.

## Installation

You can install the package by copying the command below into `R`.
``` r
remotes::install_github(repo="phillipbvetter/ctsmrTMB", dependencies=TRUE)
```

Note that `ctsmrTMB` depends on the `TMB` package. Windows users must have Rtools intalled, and Mac users must have command line tools for C++ compilers etc. You can find the GitHub for TMB [here](https://github.com/kaskr/adcomp) and installation instructions [here](https://github.com/kaskr/adcomp/wiki/Download)

## How to get started
You can visit the package [webpage](https://phillipbvetter.github.io/ctsmrTMB/index.html) and browse the vignettes for example uses, in particular see [Getting Started](https://phillipbvetter.github.io/ctsmrTMB/articles/ctsmrTMB.html).

## Help
You can access the documentation for the available methods via
``` r
?ctsmrTMB
```
This methods documentation is also available on the [homepage](https://phillipbvetter.github.io/ctsmrTMB/reference/ctsmrTMB.html).

## Example Usage

```r
library(ctsmrTMB)
library(ggplot2)
library(patchwork)

############################################################
# Data simulation
############################################################

# Simulate data using Euler Maruyama
set.seed(10)
pars = c(theta=10, mu=1, sigma_x=1, sigma_y=1e-2)
# 
dt.sim = 1e-3
t.sim = seq(0,1,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
x = 3
for(i in 1:(length(t.sim)-1)) {
  x[i+1] = x[i] + pars[1]*(pars[2]-x[i])*dt.sim + pars[3]*dw[i]
}

# Extract observations and add noise
dt.obs = 1e-2
t.obs = seq(0,1,by=dt.obs)
y = x[t.sim %in% t.obs] + pars[4] * rnorm(length(t.obs))

# Create data
.data = data.frame(
  t = t.obs,
  y = y
)

############################################################
# Model creation and estimation
############################################################

# Create model object
obj = ctsmrTMB$new()

# Set name of model (and the created .cpp file)
obj$set_modelname("ornstein_uhlenbeck")

# Set path where generated C++ files are saved.
# This will create a cppfiles folder in your current working directory if it doesnt exist
obj$set_cppfile_directory("cppfiles")

# Add system equations
obj$add_systems(
  dx ~ theta * (mu-x) * dt + sigma_x*dw
)

# Add observation equations
obj$add_observations(
  y ~ x
)

# Set observation equation variances
obj$add_observation_variances(
  y ~ sigma_y^2
)

# Specify algebraic relations
obj$add_algebraics(
  theta   ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)

# Specify parameter initial values and lower/upper bounds in estimation
obj$add_parameters(
  logtheta    = log(c(init = 1, lower=1e-5, upper=50)),
  mu          = c(init=1.5, lower=0, upper=5),
  logsigma_x  = log(c(init= 1e-1, lower=1e-10, upper=10)),
  logsigma_y  = log(c(init=1e-1, lower=1e-10, upper=10))
)

# Set initial state mean and covariance
obj$set_initial_state(x[1], 1e-1*diag(1))

# Carry out estimation using extended kalman filter method with stats::nlminb as optimizer
fit <- obj$estimate(data=.data, method="ekf", ode.solver="rk4", use.hessian=TRUE)

# Check parameter estimates against truth
cbind(c(exp(x[1]),x[2],exp(x[3]),exp(x[4])), pars)

# plot one-step predictions, simulated states and observations
t.est = fit$states$mean$prior$t
x.mean = fit$states$mean$prior$x
x.sd = fit$states$sd$prior$x
ggplot() +
  geom_ribbon(aes(x=t.est, ymin=x.mean-2*x.sd, ymax=x.mean+2*x.sd),fill="grey", alpha=0.9) +
  geom_line(aes(x=t.est, x.mean),col="steelblue",lwd=1) +
  geom_line(aes(x=t.sim,y=x)) + 
  geom_point(aes(x=t.obs,y=y),col="tomato",size=2) +
  theme_minimal()


# Plot one-step-ahead residual analysis
plot(fit, use.ggplot=T)

# Predict to obtain k-step-ahead predictions to see model forecasting ability
pred = obj$predict(data=.data, k.ahead=10, method="ekf", ode.solver="rk4")

# Plot all 10-step predictions against data
pred10step = pred[pred$k.ahead==10,]
ggplot() +
  geom_point(aes(x=t.obs,y=y),color="steelblue") +
  geom_point(aes(x=pred10step$t.j,pred10step$x),color="tomato") +
  geom_errorbar(aes(x=pred10step$t.j, 
                    ymin=pred10step$x-2*sqrt(pred10step$var.x), 
                    ymax=pred10step$x+2*sqrt(pred10step$var.x)),color="tomato") +
  theme_minimal()
```


