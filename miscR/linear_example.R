# We consider the linear ornstein-uhlenbeck
# dx = theta * ( mu - x ) * dt + sigmaX * dw
# y = x + e ~ N(0,sigmaY)

# Initialize
rm(list=ls())

# Libraries
library(devtools)
library(Deriv)
library(TMB)

# Set directory
this_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(this_path))

# Loading functions package style 
pack = as.package("../R")
load_all(pack)

# Loading functions directly
# source("../R/makeNLL.R")
# source("../R/write_linearExact_cpp.R")
# source("../R/write_TMB_cpp.R")
# source("../R/write_ExtendedKalman_cpp.R")
# source("../R/IsSystemLinear.R")
# source("../R/hat2pow.R")

# Construct model list
model = list()

model$modelname = "linear_example"

model$sdeEq = list()
model$sdeEq[[1]] = expression( dX ~ exp(logtheta) * (mu - X) * dt + exp(logsigmaX) * dw1 )
# Note: Currently not allowed to use dw, must be numbered e.g. dw1, dw2,...

model$obsEq = list()
model$obsEq[[1]] = expression(Y ~ X)

model$obsVar = list()
model$obsVar[[1]] = expression(exp(logsigmaY))

# Construct data list (different for kalman filter and for TMB)

# Simulation:
# parameters
theta  = 1
mu     = 0.5
sigmaX = 1e-1
sigmaY = 1e-2

# initial state
X0 = 1

# simulation functions
f     = function(x) theta * (mu - x)
g     = function(x) sigmaX
euler = function(x, f, g, tvec, dB=NULL){
  X <- numeric(length(tvec))
  X[1] <- X0
  dt <- diff(tvec)
  if(is.null(dB)) dB <- rnorm(length(dt), sd=sqrt(dt))
  for(i in 1:(length(tvec)-1))
    X[i+1] <- X[i] + f(X[i])*dt[i] + g(X[i])*dB[i]
  return(X)
}

# time variables
Tstep  = 0.1     # Time step for Euler scheme
Tfinal = 50      # Duration of simulation
Tobs   = 1       # Sample time interval. Should be divisible with Tsim

# Simulate trajectory
tsim <- seq(0, Tfinal, Tstep)
Xsim <- euler(X0, f, g, tsim)

# Measurements with sample interval Tobs
iobs <- seq(1, length(tsim), round(Tobs/Tstep))

# generate random measurements
Y.tmb = rnorm(length(iobs), mean=(Xsim[iobs]), sd = sigmaY)
Y.kalman = rep(NA,length(tsim))
Y.kalman[iobs] = Y.tmb
# Note: currently TMB-style uses iobs, and Kalman Filter checks for NA-values instead such that
# it is not necessary to provide iobs. Let's discuss these two approaches.

plot(tsim,Xsim,type="l")
points(tsim[iobs],Y.tmb,type="p",pch=16,col="red")

################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################

data.tmb = list()
data.tmb$t = tsim
data.tmb$Y = Y.tmb
data.tmb$iobsY = iobs-1
data.tmb$X = Xsim
data.tmb$pars = list(
  logtheta = log(0.1),
  mu = 2,
  logsigmaX = log(1),
  logsigmaY = log(1)
)

################ KALMAN DATA ################
################ KALMAN DATA ################
################ KALMAN DATA ################

data.kalman = list()
data.kalman$t = tsim
data.kalman$Y = Y.kalman
data.kalman$X0 = as.vector(0.5) #initial state expectation
data.kalman$P0 = as.matrix(0.1) #initial state covariance
data.kalman$dt = diff(tsim)[1]/10 #time-step for integration in ODE for 1st and 2nd order moments
data.kalman$pars = list(
  logtheta = log(0.1),
  mu = 2,
  logsigmaX = log(1),
  logsigmaY = log(1)
)

################################################################################
# CASE 1 - RESULTS USING TMB-STYLE
################################################################################

# construct nll for tmb-style
nll.tmb = makeNLL(model,data.tmb,method="TMB")

# Estimate parameters and latent variables
time.tmb = system.time(opt.tmb <- nlminb(nll.tmb$par,nll.tmb$fn,nll.tmb$gr))
sdr <- sdreport(nll.tmb)

# Get predictions of states with std.dev.
pred <- summary(sdr,"random")
Xpred <- pred[,1]
Xsd <- pred[,2]

# Setup plot
plot(tsim, Xsim, type="n", xlab="Time t", ylab="State x")

# Confidence region for states
polygon(c(tsim, rev(tsim)), c(Xpred+1.96*Xsd,rev(Xpred-1.96*Xsd)), col="grey", border=NA)

# Plot true path
lines(tsim, Xsim, type="l")

# Add predictions
lines(tsim,Xpred,pch=16,col="red")

# Add measurements
points(tsim[iobs],Y.tmb)

legend("topright",legend=c("True","Measured","Smoothed"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# CASE 2 - RESULTS USING KALMAN FILTER-STYLE
################################################################################

# so cpp files do not overwrite
model.kalman = model
model.kalman$modelname = "linear_example_kalman"
# construct nll for kalman filter-style
nll.kalman = makeNLL(model.kalman,data.kalman,method="kalman")

# Estimate parameters and latent variables
time.kalman = system.time(opt.kalman <- nlminb(nll.kalman$par,nll.kalman$fn,nll.kalman$gr))

Xpred = unlist( nll.kalman$report()$xPost )
Xsd   = sqrt(unlist( nll.kalman$report()$pPost ))

# Setup plot
plot(tsim, Xsim, type="n", xlab="Time t", ylab="State x")

# Confidence region for states
polygon(c(tsim, rev(tsim)), c(Xpred+1.96*Xsd,rev(Xpred-1.96*Xsd)), col="grey", border=NA)

# Plot true path
lines(tsim, Xsim, type="l")

# Add predictions
lines(tsim,Xpred,pch=16,col="red")

# Add measurements
points(tsim[iobs],Y.tmb)

legend("topright",legend=c("True","Measured","Posterior"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# CASE 3 - RESULTS USING EXACT-STYLE
################################################################################

# so cpp files do not overwrite
model.exact = model
model.exact$modelname = "linear_example_exact"
# construct nll for exact-style
nll.exact = makeNLL(model.exact,data.tmb,exact=TRUE)
# Estimate parameters and latent variables
time.exact = system.time(opt.exact <- nlminb(nll.exact$par,nll.exact$fn,nll.exact$gr))
sdr <- sdreport(nll.exact)

# Get predictions of states with std.dev.
pred <- summary(sdr,"random")
Xpred <- pred[,1]
Xsd <- pred[,2]

# Setup plot
plot(tsim, Xsim, type="n", xlab="Time t", ylab="State x")

# Confidence region for states
polygon(c(tsim, rev(tsim)), c(Xpred+1.96*Xsd,rev(Xpred-1.96*Xsd)), col="grey", border=NA)

# Plot true path
lines(tsim, Xsim, type="l")

# Add predictions
lines(tsim,Xpred,pch=16,col="red")

# Add measurements
points(tsim[iobs],Y.tmb)

legend("topright",legend=c("True","Measured","Smoothed"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# COMPARE RESULTS
################################################################################

# Run-time for optimization
print(times <- rbind(time.tmb, 
                     time.kalman,
                     time.exact)
      )

# Parameter values
log2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4]))
print(pars <- rbind(True = c(theta=1,mu=1,sigmaX=1e-1,sigmaY=1e-2),
                    TMB = log2real(opt.tmb$par),
                    Kalman = log2real(opt.kalman$par),
                    Exact = log2real(opt.exact$par))
      )

# Objective value at optimum
print(nll <- rbind(TMB = opt.tmb$objective,
                   Kalman = opt.kalman$objective,
                   Exact = opt.exact$objective)
      )
