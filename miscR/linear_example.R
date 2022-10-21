# We consider the linear ornstein-uhlenbeck
# dx = theta * ( mu - x ) * dt + sigmaX * dw
# y = x + e ~ N(0,sigmaY)

# Initialize
rm(list=ls())

# Libraries
library(devtools)
library(Deriv)
library(TMB)
# library(ctsmr)


this_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(this_path))

# Loading functions package style 
pack = as.package("../R")
load_all(pack)

# Construct model list
model = list()

model$sdeEq = list()
model$sdeEq[[1]] = expression( dx ~ theta * (mu - x) * dt + sigma_x * dw )
# Note: Currently not allowed to use dw, must be numbered e.g. dw1, dw2,...

model$obsEq = list()
model$obsEq[[1]] = expression(y ~ x)

model$obsVar = list()
model$obsVar[[1]] = expression(sigma_y)

model$algeqs = list()
model$algeqs[["theta"]] = expression(exp(logtheta))
model$algeqs[["sigma_x"]] = expression(exp(logsigma_x))
model$algeqs[["sigma_y"]] = expression(exp(logsigma_y))

# Construct data list (different for kalman filter and for TMB)

# Simulation:
# parameters
theta  = 1
mu     = 1
sigmaX = 1e-1   
sigmaY = 1e-2^2 #obs variance

# initial state
X0 = 1.5

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
Y.data = rep(NA,length(tsim))
Y.data[iobs] = rnorm(length(iobs), mean=(Xsim[iobs]), sd = sqrt(sigmaY))

plot(tsim,Xsim,type="l")
points(tsim[iobs],Y.data[iobs],type="p",pch=16,col="red")

################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################

theta_init = 0.1
mu_init = mean(Y.data,na.rm=T)
sigma_x_init = 0.5
sigma_y_init = var(Y.data,na.rm=T)

data.tmb = list()
data.tmb$t = tsim
data.tmb$y = Y.data
data.tmb$x = rep(mean(Y.data,na.rm=T),length(Xsim))
data.tmb$pars = list(
  logtheta = log(theta_init),
  mu = mu_init,
  logsigma_x = log(sigma_x_init),
  logsigma_y = log(sigma_y_init)
)

################ KALMAN DATA ################
################ KALMAN DATA ################
################ KALMAN DATA ################

data.kalman = list()
data.kalman$t = tsim
data.kalman$y = Y.data
data.kalman$X0 = as.vector(Y.data[1]) #initial state expectation
data.kalman$P0 = as.matrix(0.1^2) #initial state covariance
data.kalman$dt = diff(tsim)[1]/10 #time-step for integration in ODE for 1st and 2nd order moments
data.kalman$pars = list(
  logtheta = log(theta_init),
  mu = mu_init,
  logsigma_x = log(sigma_x_init),
  logsigma_y = log(sigma_y_init)
)

################################################################################
# RESULTS USING TMB
################################################################################

model$modelname = "nll_tmb"

# construct nll for tmb-style
nll.tmb = makeNLL(model,data.tmb,method="tmb",compile=F,silent=T)

# Estimate parameters and latent variables
time.tmb = mark(
  opt.tmb <- nlminb(nll.tmb$par,nll.tmb$fn,nll.tmb$gr),iterations=10
)$median
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
points(tsim[iobs],Y.data[iobs])

legend("topright",legend=c("True","Measured","Smoothed"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# RESULTS USING TMB EXACT
################################################################################

# so cpp files do not overwrite
model.exact = model
model.exact$modelname = "nll_tmbexact"
# construct nll for exact-style
nll.exact = makeNLL(model.exact,data.tmb,method="tmb_exact",compile=F,,silent=T)
# Estimate parameters and latent variables
time.exact = mark(
  opt.exact <- nlminb(nll.exact$par,nll.exact$fn,nll.exact$gr),iterations=10
)$median
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
points(tsim[iobs],Y.data[iobs])

legend("topright",legend=c("True","Measured","Smoothed"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# CASE 2 - RESULTS USING EXTENDED KALMAN FILTER
################################################################################

# so cpp files do not overwrite
model.ekf = model
model.ekf$modelname = "nll_ekf"
# construct nll for kalman filter-style
nll.ekf = makeNLL(model.ekf,data.kalman,method="ekf",compile=F,,silent=T)

# Estimate parameters and latent variables
time.ekf = mark(
  opt.ekf <- nlminb(nll.ekf$par,nll.ekf$fn,nll.ekf$gr),iterations=10
)$median
# time.ekf = system.time(opt.ekf <- nlminb(nll.ekf$par,nll.ekf$fn,nll.ekf$gr,nll.ekf$he),gcFirst=T)

Xpred = unlist( nll.ekf$report()$xPost )
Xsd   = sqrt(unlist( nll.ekf$report()$pPost ))

# Setup plot
plot(tsim, Xsim, type="n", xlab="Time t", ylab="State x")

# Confidence region for states
polygon(c(tsim, rev(tsim)), c(Xpred+1.96*Xsd,rev(Xpred-1.96*Xsd)), col="grey", border=NA)

# Plot true path
lines(tsim, Xsim, type="l")

# Add predictions
lines(tsim,Xpred,pch=16,col="red")

# Add measurements
points(tsim[iobs],Y.data[iobs])

legend("topright",legend=c("True","Measured","Posterior"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# CASE 2 - RESULTS USING UNSCENTED KALMAN FILTER
################################################################################

# so cpp files do not overwrite
model.ukf = model
model.ukf$modelname = "nll_ukf"
# construct nll for kalman filter-style

nll.ukf = makeNLL(model.ukf,data.kalman,method="ukf",compile=F,silent=T)

# Estimate parameters and latent variables
time.ukf = mark(
  opt.ukf <- nlminb(nll.ukf$par,nll.ukf$fn,nll.ukf$gr),iterations=10
)$median
# time.ukf = system.time(opt.ukf <- nlminb(nll.ukf$par,nll.ukf$fn,nll.ukf$gr,nll.ukf$he),gcFirst=T)

Xpred = unlist( nll.ukf$report()$xPost )
Xsd   = sqrt(unlist( nll.ukf$report()$pPost ))

# Setup plot
plot(tsim, Xsim, type="n", xlab="Time t", ylab="State x")

# Confidence region for states
polygon(c(tsim, rev(tsim)), c(Xpred+1.96*Xsd,rev(Xpred-1.96*Xsd)), col="grey", border=NA)

# Plot true path
lines(tsim, Xsim, type="l")

# Add predictions
lines(tsim,Xpred,pch=16,col="red")

# Add measurements
points(tsim[iobs],Y.data[iobs])

legend("topright",legend=c("True","Measured","Posterior"),lty=c("solid",NA,"solid"),pch=c(NA,1,NA),col=c("black","black","red"))

################################################################################
# USING CTSM-R
################################################################################

# # Empty model object
# obj = ctsm$new()
# 
# # Add a system equation
# obj$addSystem(dx ~ exp(logtheta)*(mu-x)*dt + exp(logsigma)*dw)
# 
# # Add an observation equation
# obj$addObs(y ~ x)
# 
# # Set the observation variance
# obj$setVariance(y ~ exp(logsigmaY))
# 
# # Set initial model parameters
# obj$setParameter(logtheta  = c(init=log(0.1),lb=log(1e-6),  ub=log(10)))
# obj$setParameter(mu        = c(init=2,       lb=-10,        ub=10))
# obj$setParameter(logsigma  = c(init=log(1),  lb=log(1e-6),  ub=log(2)))
# obj$setParameter(logsigmaY = c(init=log(1),  lb=log(1e-6),  ub=log(2)))
# 
# # Set initial state values
# obj$setParameter(x = c(init=1.5,lb=0.5,ub=2.5))
# 
# dat = data.frame(t=tsim,y=Y.kalman)
# 
# time.ctsmr = system.time(fit <- obj$estimate(dat),gcFirst=T)

################################################################################
# COMPARE RESULTS
################################################################################

# # Run-time for optimization
# print(times <- rbind(TMB = time.tmb, 
#                      Kalman = time.kalman,
#                      Exact = time.exact,
#                      CTSMR = time.ctsmr)[,c(1,3)]
#       )
# 
# # Parameter values
# log2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4]))
# print(pars <- rbind(True = c(theta,mu,sigmaX,sigmaY),
#                     TMB = log2real(opt.tmb$par),
#                     Kalman = log2real(opt.kalman$par),
#                     Exact = log2real(opt.exact$par),
#                     CTSMR = log2real(fit$xm[c("logtheta","mu","logsigma","logsigmaY")]))
#       )
# 
# # Objective value at optimum
# print(nll <- rbind(TMB = opt.tmb$objective,
#                    Kalman = opt.kalman$objective,
#                    Exact = opt.exact$objective,
#                    CTSMR = fit$f)
#       )


################################################################################
# COMPARE RESULTS WITHOUT CTSM-R
################################################################################

# Run-time for optimization
print(times <- rbind(tmb = time.tmb, 
                     tmb_exact = time.exact,
                     ekf = time.ekf,
                     ukf = time.ukf)
)

# Parameter values
log2real = function(x) c(exp(x[1]),x[2],exp(x[3]),exp(x[4]))
print(pars <- rbind(true = c(theta,mu,sigmaX,sigmaY),
                    tmb = log2real(opt.tmb$par),
                    tmb_exact = log2real(opt.exact$par),
                    ekf = log2real(opt.ekf$par),
                    ukf = log2real(opt.ukf$par))
)

# Objective value at optimum
print(nll <- rbind(tmb = opt.tmb$objective,
                   tmb_exact = opt.exact$objective,
                   ekf = opt.ekf$objective,
                   ukf = opt.ukf$objective)
)
