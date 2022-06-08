# OBS: You can view results at the bottom of this file instead of running the code
# by loading the saved optdata file!

# We are considering the linear ornstein-uhlenbeck process with gaussian observation noise
# dx = theta * ( mu - x ) * dt + sigmaX * dw
# y = x + e, where e ~ N(0,sigmaY)

# Initialize
rm(list=ls())

# Libraries
library(devtools)
library(Deriv)
library(TMB)
library(ctsmr)

# Set directory
this_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(this_path))

# Loading functions package style 
pack = as.package("../R")
load_all(pack)

# Construct model list
model = list()

model$sdeEq = list()
model$sdeEq[[1]] = expression( dX ~ theta * (mu - X) * dt + sigmaX * dw1 )
# Note: Currently not allowed to use dw, must be numbered e.g. dw1, dw2,...

model$obsEq = list()
model$obsEq[[1]] = expression(Y ~ X)

model$obsVar = list()
model$obsVar[[1]] = expression(sigmaY)

model$algeqs = list()
model$algeqs[["theta"]] = expression(exp(logtheta))
model$algeqs[["sigmaX"]] = expression(exp(logsigmaX))
model$algeqs[["sigmaY"]] = expression(exp(logsigmaY))

# Simulation:
# parameters
theta  = 1
mu     = 0.5
sigmaX = 1e-1
sigmaY = 1e-2

# initial guesses
theta.init  = 0.1
mu.init     = 2
sigmaX.init = 1
sigmaY.init = 1

# initial state
X0 = 0.6

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

# Time vector
tsim <- seq(0, Tfinal, Tstep)

# Measurements with sample interval Tobs
iobs <- seq(1, length(tsim), round(Tobs/Tstep))

# Simulate trajectories
N = 100
Xsim = matrix(nrow=N,ncol=length(tsim))
Y.tmb = matrix(nrow=N,ncol=length(iobs))
Y.kalman = Xsim
for(i in 1:N) {
  Xsim[i,] <- euler(X0, f, g, tsim)
  Y.tmb[i,] <- rnorm(length(iobs), mean=(Xsim[i,iobs]), sd = sigmaY)
  Y.kalman[i,] <- rep(NA,length(tsim))
  Y.kalman[i,iobs] = Y.tmb[i,]
  
}

################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################

data.tmb = list()
data.tmb$t = tsim
data.tmb$iobsY = iobs-1
data.tmb$Y = Y.tmb[1,]
data.tmb$X = rep(mean(data.tmb$Y),length(data.tmb$Y))

pars.tmb = list()
pars.tmb$logtheta = log(theta.init)
pars.tmb$mu = mu.init
pars.tmb$logsigmaX = log(sigmaX.init)
pars.tmb$logsigmaY = log(sigmaY.init)

data.tmb$pars = pars.tmb
################ KALMAN DATA ################
################ KALMAN DATA ################
################ KALMAN DATA ################

data.kalman = list()
data.kalman$t = tsim
data.kalman$dt = diff(tsim)[1]/10 #time-step for integration in ODE for 1st and 2nd order moments
data.kalman$Y = Y.kalman[1,]
data.kalman$X0 = mean(data.kalman$Y,na.rm=T)
data.kalman$P0 = sd(data.kalman$Y,na.rm=T)

pars.kalman = list()
pars.kalman$logtheta = log(theta.init)
pars.kalman$mu = mu.init
pars.kalman$logsigmaX = log(sigmaX.init)
pars.kalman$logsigmaY = log(sigmaY.init)

data.kalman$pars = pars.kalman

# optimization upper and lower bounds
lb = c( log(1e-5), -10, log(1e-5), log(1e-5) )
ub = c( log(10), 10, log(2), log(2) )

################################################################################
# CASE 1 - RESULTS USING TMB-STYLE
################################################################################

model$modelname = "linear_example_TMB"
if(!file.exists("linear_example_TMB.cpp")){
  write_TMB_cpp(model,data.tmb)
  compile("linear_example_TMB.cpp")
}
dyn.load(dynlib("linear_example_TMB"))

tmb.opt.results = list()
tmb.pars.results = list()
tmb.cor.results = list()
tmb.time.results = list()


for(i in 1:N){
  data.tmb$Y = Y.tmb[i,]
  pars.tmb$X = rep(mean(data.tmb$Y),length(data.tmb$t))
  # Return neg. log likelihood
  nll = MakeADFun(data.tmb, pars.tmb, random="X", DLL="linear_example_TMB",silent=T)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
    tmb.time.results[[i]] = system.time(opt <- nlminb(nll$par,nll$fn,nll$gr,lower=lb,upper=ub))
    sdr <- sdreport(nll)
    tmb.opt.results[[i]] = opt
    tmb.pars.results[[i]] = sdr$par.fixed
    tmb.cor.results[[i]] = cov2cor(sdr$cov.fixed)
    },
    error = function(e) print(e),
    warning = function(w) print(w)
    )
  print(i)
}

################################################################################
# CASE 2 - RESULTS USING EXACT-STYLE
################################################################################

# so cpp files do not overwrite
model$modelname = "linear_example_exact"
if(!file.exists("linear_example_exact.cpp")){
  write_linearExact_cpp(model,data.tmb)
  compile("linear_example_exact.cpp")
}
dyn.load(dynlib("linear_example_exact"))

exact.opt.results = list()
exact.pars.results = list()
exact.cor.results = list()
exact.time.results = list()

for(i in 1:N){
  data.tmb$Y = Y.tmb[i,]
  pars.tmb$X = rep(mean(data.tmb$Y),length(data.tmb$t))
  # Return neg. log likelihood
  nll = MakeADFun(data.tmb, pars.tmb, random="X", DLL="linear_example_exact",silent=T)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
    exact.time.results[[i]] = system.time(opt <- nlminb(nll$par,nll$fn,nll$gr,lower=lb,upper=ub))
    sdr <- sdreport(nll)
    exact.opt.results[[i]] = opt
    exact.pars.results[[i]] = opt$par
    exact.cor.results[[i]] = cov2cor(sdr$cov.fixed)
    },
    error = function(e) print(e),
    warning = function(w) print(w)
  )
  print(i)
}


################################################################################
# CASE 3 - RESULTS USING KALMAN FILTER-STYLE
################################################################################

model$modelname = "linear_example_kalman"
if(!file.exists("linear_example_kalman.cpp")){
  write_ExtendedKalman_cpp(model,data.kalman)
  compile("linear_example_kalman.cpp")
}
dyn.load(dynlib("linear_example_kalman"))

kalman.pars.results = list()
kalman.cor.results = list()
kalman.opt.results = list()
kalman.time.results = list()

for(i in 1:N){
  data.kalman$Y = Y.kalman[i,]
  data.kalman$X0 = mean(data.kalman$Y,na.rm=T)
  data.kalman$P0 = matrix(sd(data.kalman$Y,na.rm=T))
  # Return neg. log likelihood
  nll = MakeADFun(data.kalman, pars.kalman, DLL="linear_example_kalman",silent=T)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
      kalman.time.results[[i]] = system.time(opt <- nlminb(nll$par,nll$fn,nll$gr,nll$he,lower=lb,upper=ub))
      kalman.opt.results[[i]] = opt
      kalman.pars.results[[i]] = opt$par
      kalman.cor.results[[i]] = cov2cor(solve(nll$he(opt$par)))
    },
    error = function(e) print(e),
    warning = function(w) {print(w);print(solve(nll$he(opt$par)))}
  )
  print(i)
}

################################################################################
# CASE 4 - USING CTSM-R
################################################################################

# Empty model object
obj = ctsm$new()
# Add a system equation
obj$addSystem(dx ~ exp(logtheta)*(mu-x)*dt + exp(logsigmaX)*dw)
# Add an observation equation
obj$addObs(y ~ x)
# Set the observation variance
obj$setVariance(y ~ exp(logsigmaY))
# Set initial model parameters
obj$setParameter(logtheta  = c(init=log(theta.init),  lb=-5,    ub=ub[1]))
obj$setParameter(mu        = c(init=mu.init,          lb=lb[2], ub=ub[2]))
obj$setParameter(logsigmaX  = c(init=log(sigmaX.init),lb=-5,    ub=ub[3]))
obj$setParameter(logsigmaY = c(init=log(sigmaY.init), lb=-5,    ub=ub[4]))

ctsm.opt.results = list()
ctsm.pars.results = list()
ctsm.cor.results = list()
ctsm.time.results = list()

for(i in 1:N){
  data.ctsm = data.frame(t=tsim,y=Y.kalman[i,])
  x.init = mean(data.ctsm$y,na.rm=T) + c(0,-2,2)*sd(data.ctsm$y,na.rm=T)
  obj$setParameter(x = c(init=x.init[1],lb=x.init[2],ub=x.init[3]))
  out = tryCatch(
    {
      ctsm.time.results[[i]] = system.time(fit <- obj$estimate(data.ctsm))
      ctsm.opt.results[[i]] = fit
      ctsm.pars.results[[i]] = fit$xm[-1][c("logtheta","mu","logsigmaX","logsigmaY")]
      ctsm.cor.results[[i]] = fit$corr[2:5,2:5]
    },
    error = function(e) print(e),
    warning = function(w) print(w)
  )
  print(i)
}

################################################################################
# COMPARE RESULTS
################################################################################

# optdata = list(tmb.opt.results=tmb.opt.results, tmb.pars.results=tmb.pars.results, tmb.cor.results=tmb.cor.results,
#                tmb.time.results=tmb.time.results, exact.opt.results=exact.opt.results, exact.pars.results=exact.pars.results,
#                exact.cor.results=exact.cor.results, exact.time.results=exact.time.results,kalman.opt.results=kalman.opt.results,
#                kalman.pars.results=kalman.pars.results, kalman.cor.results=kalman.cor.results, kalman.time.results=kalman.time.results,
#                ctsm.opt.results=ctsm.opt.results, ctsm.pars.results=ctsm.pars.results, ctsm.cor.results=ctsm.cor.results,
#                ctsm.time.results=ctsm.time.results)
# save(optdata,file="optdata.Rda")

if(!(exists("tmb.pars.results") & exists("exact.pars.results") & exists("kalman.pars.results") & exists("ctsm.pars.results"))){
  load("optdata.Rda")
  for(name in names(optdata)){
    assign(name,optdata[[name]])
  }
}

# true parameter values
theta  = 1
mu     = 0.5
sigmaX = 1e-1
sigmaY = 1e-2

# check convergence status
sapply(tmb.opt.results,function(x) x$message)
sapply(exact.opt.results,function(x) x$message)
sapply(kalman.opt.results,function(x) x$message)
sapply(ctsm.opt.results,function(x) x$message)

# boxplot for theta
theta.data=
  cbind(
  tmb=sapply(tmb.pars.results,function(x) exp(x[1])),
  exact=sapply(exact.pars.results,function(x) exp(x[1])),
  kalman=sapply(kalman.pars.results,function(x) exp(x[1])),
  ctsm=sapply(ctsm.pars.results,function(x) exp(x[1]))
)
boxplot(theta.data)
abline(h=theta,lw=2,col="red")

# boxplot for mu
mu.data=
  cbind(
    tmb=sapply(tmb.pars.results,function(x) x[2]),
    exact=sapply(exact.pars.results,function(x) x[2]),
    kalman=sapply(kalman.pars.results,function(x) x[2]),
    ctsm=sapply(ctsm.pars.results,function(x) x[2])
  )
boxplot(mu.data)
abline(h=mu,lw=2,col="red")

# boxplot for sigmaX
sigmaX.data=
  cbind(
    tmb=sapply(tmb.pars.results,function(x) exp(x[3])),
    exact=sapply(exact.pars.results,function(x) exp(x[3])),
    kalman=sapply(kalman.pars.results,function(x) exp(x[3])),
    ctsm=sapply(ctsm.pars.results,function(x) exp(x[3]))
  )
boxplot(sigmaX.data)
abline(h=sigmaX,lw=2,col="red")

# boxplot for sigmaY
sigmaY.data=
  cbind(
    tmb=sapply(tmb.pars.results,function(x) exp(x[4])),
    exact=sapply(exact.pars.results,function(x) exp(x[4])),
    kalman=sapply(kalman.pars.results,function(x) exp(x[4])),
    ctsm=sapply(ctsm.pars.results,function(x) exp(x[4]))
  )
boxplot(sigmaY.data)
abline(h=sigmaY,lw=2,col="red")

time.comparison = cbind(
  tmb = sapply(tmb.time.results,function(x) x[1]),
  exact = sapply(exact.time.results,function(x) x[1]),
  kalman = sapply(kalman.time.results,function(x) x[1]),
  ctsm = sapply(ctsm.time.results,function(x) x[1])
)
boxplot(time.comparison)
