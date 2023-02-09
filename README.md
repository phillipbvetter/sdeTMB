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

## Example Usage

``` r
# Create model object
obj = sdemTMB$new()

# Modelname
obj$set_modelname("ekf_model")

# Method (choose "ekf", "ukf" or "tmb")
obj$set_method("ekf")

# Set location path for C++ file generation
# obj$set_path("~/")

# System
obj$add_systems(
  dx ~ theta * (mu-x) * dt + sigma_x*dw + sigma_x*dw2,
  du ~ sigma_u * dw
)
# Obs 
obj$add_observations(
  y ~ x,
  w ~ u
)
# Obs Variance
obj$add_observation_variances(
  y ~ sigma_y^2,
  w ~ sigma_w^2
)
# Algebra
obj$add_algebraics(
  theta ~ exp(logtheta),
  mu ~ mu0 + a*Qf + b*Sf + c*Qr,
  sigma_x ~ exp(logsigma_x),
  sigma_u ~ exp(logsigma_u),
  sigma_y ~ exp(logsigma_y),
  sigma_w ~ exp(logsigma_w)
)
# Inputs
obj$add_inputs(Qf,Sf,Qr)

# Parameters
obj$add_parameters(
  logtheta ~ log(c(1,1e-5,2)),
  logsigma_x ~ log(c(1e-1,1e-20,2)),
  logsigma_u ~ log(c(1e-1,1e-20,2)),
  logsigma_y ~ log(c(1e-1,1e-20,2)),
  logsigma_w ~ log(c(1e-1,1e-20,2)),
  mu0 ~ c(1,-10,10),
  a ~ c(0,-100,100),
  b ~ c(0,-100,100),
  c ~ c(0,-100,100)
)

# obj$add_constants(R~1,V~2.5)

# Initial State
# obj$set_initial_state(1,diag(1))
obj$set_initial_state(c(1,1),1e-2*diag(2))

# obj$set_lamperti("log",c("z","u"))
# obj$set_lamperti("log")
# obj$set_lamperti("log",c("z"))

set.seed(10)
t = seq(0,10,length.out=n)
data = data.frame(
  t = t ,
  Qf = rep(1,n),
  Sf = rep(1,n),
  Qr = rep(1,n),
  y = rnorm(n,mean=1,sd=0.2),
  w = rnorm(n,mean=1,sd=0.05)
)
# id1 = sample(2:(n-1),size=round(n*5/10),replace=FALSE)
# id2 = sample(2:(n-1),size=round(n*5/10),replace=FALSE)
# data$y[id1] = NA
# data$w[id2] = NA
# n1 = round(n/4)
# n2 = round(n*3/4)
# data$y[n1:n2] = NA
# data$w[n1:n2] = NA

# obj$set_timestep(0.25)
# obj$set_silence(FALSE)
# obj$use_hessian(T)

# self = obj

# return.fit = TRUE
# return.nll = FALSE
# use.hessian=FALSE
# ode.timestep = 0.25
# compile = FALSE
# silent = TRUE

fit <- obj$estimate(data, silence=T,
                    compile=FALSE)

# plot(fit, use.ggplot=F, plot.obs=2)

pred = predict(fit,
               n.step.ahead=10,
               use.simulation=TRUE,
               return.state.dispersion=TRUE,
               covariance=TRUE,
               sim.timestep=0.05,
               n.sim = 1e3,
               give.only.n.step.ahead=F)
```


