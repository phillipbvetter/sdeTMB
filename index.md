### Stochastic Differential Equations using Template Model Builder

___

Welcome to the homepage for `sdeTMB`.

`sdeTMB` is an `R` package for parameter estimation, state filtration and forecasting in stochastic state space system. 

The package is maintained by a group of researchers at the Technical University of Denmark (DTU), Department of Applied Mathematics and Computer Science, at the section for Dynamical Systems.

The package is the successor of the `ctsmr` package whose FORTRAN libraries became increasingly difficult to maintain at the section. We now instead rely on the `TMB` package which gives easy access to a series of fast and modern `C++` libraries, including automatic differentiation through `CppAD`.

The package implements a `R6` class object through which the user specifies a state space model which consists of a system of stochastic differential equations (SDEs) and a series of algebraic observation equations (AOEs). The package uses the defined model system and a user-supplied dataset of observations to write a `C++` function in the form accepted by `TMB` which calculates the negative log-likelihood. 

### Installation

---

You can install the latest version using

``` r
remotes::install_github(repo="phillipbvetter/sdeTMB", dependencies=TRUE)
```
You should also follow the installation instructions for `TMB` [here](https://github.com/kaskr/adcomp/wiki/Download) to get everything working. Windows users must have Rtools intalled, and Mac users must have Command Line Tools for various compilers. You can read more on the `TMB` GitHub [here](https://github.com/kaskr/adcomp).

### Getting Started

---

If you are a new user we recommend that you check out the [Getting Started](https://phillipbvetter.github.io/sdeTMB/articles/sdeTMB.html) vignette to start learning the simple steps of using `ctsmrTMB`.

You can find documentation for the methods [here](https://phillipbvetter.github.io/sdeTMB/reference/sdeTMB.html).

### Example Usage

---

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
x = fit$par.fixed
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










