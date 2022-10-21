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
library(bench)
# library(ctsmr)

# Set directory
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

# Simulation:
# parameters
theta  = 1
mu     = 1
sigmaX = 1e-1
sigmaY = 1e-2^2

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

# Time vector
tsim <- seq(0, Tfinal, Tstep)

# Measurements with sample interval Tobs
iobs <- seq(1, length(tsim), round(Tobs/Tstep))

# Simulate trajectories
N = 100
Xsim = matrix(nrow=N,ncol=length(tsim))
Y.data = matrix(NA,nrow=N,ncol=length(tsim))
for(i in 1:N) {
  Xsim[i,] <- euler(X0, f, g, tsim)
  Y.data[i,iobs] <- rnorm(length(iobs), mean=(Xsim[i,iobs]), sd = sqrt(sigmaY))
}

theta_init = 0.1
sigma_x_init = 1

################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################

data.tmb = list()
data.tmb$t = tsim
data.tmb$y = Y.data[1,]
data.tmb$x = rep(mean(data.tmb$y,na.rm=T),length(data.tmb$y))
data.tmb$pars = list(
  logtheta = log(0.1),
  mu = NA,
  logsigma_x = log(sigma_x_init),
  logsigma_y = NA
)

################ KALMAN DATA ################
################ KALMAN DATA ################
################ KALMAN DATA ################

data.kalman = list()
data.kalman$t = tsim
data.kalman$dt = diff(tsim)[1] #time-step for integration in ODE for 1st and 2nd order moments
data.kalman$y = Y.data[1,]
data.kalman$X0 = mean(data.kalman$y,na.rm=T)
data.kalman$P0 = matrix(var(data.kalman$y,na.rm=T))
data.kalman$pars = list(
  logtheta = log(theta_init),
  mu = NA,
  logsigma_x = log(sigma_x_init),
  logsigma_y = NA
)

# optimization upper and lower bounds
lb = c( log(1e-5), -10, log(1e-5), log(1e-5) )
ub = c( log(10), 10, log(2), log(2) )

################################################################################
# CASE 1 - RESULTS USING TMB-STYLE
################################################################################

model$modelname = "nll_tmb"

tmb.opt.results = list()
tmb.pars.results = list()
tmb.cor.results = list()
tmb.time.results = list()

for(i in 1:N){
  data.tmb$y = Y.data[i,]
  data.tmb$x = rep(mean(data.tmb$y,na.rm=T),length(data.tmb$t))
  # initial parameter guesses
  data.tmb$pars$mu = mean(data.tmb$y,na.rm=T)
  data.tmb$pars$logsigma_y = var(data.tmb$y,na.rm=T)
  # Return neg. log likelihood
  nll = makeNLL(model,data.tmb,method="tmb",compile=FALSE,silent=TRUE)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
    tmb.time.results[[i]] = mark(opt <- nlminb(nll$par,nll$fn,nll$gr,lower=lb,upper=ub),iterations=1)$median
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
model$modelname = "nll_exact"

exact.opt.results = list()
exact.pars.results = list()
exact.cor.results = list()
exact.time.results = list()

for(i in 1:N){
  data.tmb$y = Y.data[i,]
  data.tmb$x = rep(mean(data.tmb$y,na.rm=T),length(data.tmb$t))
  # initial parameter guesses
  data.tmb$pars$mu = mean(data.tmb$y,na.rm=T)
  data.tmb$pars$logsigma_y = var(data.tmb$y,na.rm=T)
  # Return neg. log likelihood
  nll = makeNLL(model,data.tmb,method="tmb_exact",compile=FALSE,silent=TRUE)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
    exact.time.results[[i]] = mark(opt <- nlminb(nll$par,nll$fn,nll$gr,lower=lb,upper=ub),iterations=1)$median
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

model$modelname = "nll_ekf"

ekf.pars.results = list()
ekf.cor.results = list()
ekf.opt.results = list()
ekf.time.results = list()

for(i in 1:N){
  data.kalman$y = Y.data[i,]
  data.kalman$X0 = mean(data.kalman$y,na.rm=T)
  data.kalman$P0 = matrix(var(data.kalman$y,na.rm=T))
  # initial parameter guesses
  data.kalman$pars$mu = mean(data.kalman$y,na.rm=T)
  data.kalman$pars$logsigma_y = var(data.kalman$y,na.rm=T)
  # Return neg. log likelihood
  nll = makeNLL(model,data.kalman,method="ekf",compile=FALSE,silent=TRUE)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
      ekf.time.results[[i]] = mark(opt <- nlminb(nll$par,nll$fn,nll$gr,nll$he,lower=lb,upper=ub),iterations=1)$median
      ekf.opt.results[[i]] = opt
      ekf.pars.results[[i]] = opt$par
      ekf.cor.results[[i]] = cov2cor(solve(nll$he(opt$par)))
    },
    error = function(e) print(e),
    warning = function(w) {print(w);print(solve(nll$he(opt$par)))}
  )
  print(i)
}

################################################################################
# CASE 3 - RESULTS USING UNSCENTED KALMAN FILTER-STYLE
################################################################################

model$modelname = "nll_ukf"

ukf.pars.results = list()
ukf.cor.results = list()
ukf.opt.results = list()
ukf.time.results = list()

for(i in 1:N){
  data.kalman$y = Y.data[i,]
  data.kalman$X0 = mean(data.kalman$y,na.rm=T)
  data.kalman$P0 = matrix(var(data.kalman$y,na.rm=T))
  # initial parameter guesses
  data.kalman$pars$mu = mean(data.kalman$y,na.rm=T)
  data.kalman$pars$logsigma_y = var(data.kalman$y,na.rm=T)
  # Return neg. log likelihood
  nll = makeNLL(model,data.kalman,method="ukf",compile=FALSE,silent=TRUE)
  # Estimate parameters and latent variables
  out <- tryCatch(
    {
      ukf.time.results[[i]] = mark(opt <- nlminb(nll$par,nll$fn,nll$gr,nll$he,lower=lb,upper=ub),iterations=1)$median
      ukf.opt.results[[i]] = opt
      ukf.pars.results[[i]] = opt$par
      ukf.cor.results[[i]] = cov2cor(solve(nll$he(opt$par)))
    },
    error = function(e) print(e),
    warning = function(w) {print(w);print(solve(nll$he(opt$par)))}
  )
  print(i)
}

################################################################################
# JULIA ADAPTIVE EKF
################################################################################

# library(JuliaCall)
# julia_setup(JULIA_HOME = "/Applications/Julia-1.8.app/Contents/Resources/julia/bin")
# 
# julia.pars.results = list()
# julia.cor.results = list()
# julia.opt.results = list()
# julia.time.results = list()
# 
# for(i in 1:N){
#   data.kalman$y = Y.data[i,]
#   data.kalman$X0 = mean(data.kalman$y,na.rm=T)
#   data.kalman$P0 = matrix(var(data.kalman$y,na.rm=T))
#   # initial parameter guesses
#   data.kalman$pars$mu = mean(data.kalman$y,na.rm=T)
#   data.kalman$pars$logsigma_y = var(data.kalman$y,na.rm=T)
#   # Return neg. log likelihood
#   julia.time = mark( nll <- makeNLL(model,data.kalman,method="kalman_adaptive",compile=FALSE,silent=TRUE),iterations=1)$median
#   # Estimate parameters and latent variables
#   julia.time.results[[i]] = julia.time
#   julia.opt.results[[i]] = nll$opt
#   julia.pars.results[[i]] = nll$opt$par
#   print(i)
# }

################################################################################
# CASE 4 - USING CTSM-R
################################################################################

# # Empty model object
# obj = ctsm$new()
# # Add a system equation
# obj$addSystem(dx ~ exp(logtheta)*(mu-x)*dt + exp(logsigmaX)*dw)
# # Add an observation equation
# obj$addObs(y ~ x)
# # Set the observation variance
# obj$setVariance(y ~ exp(logsigmaY))
# # Set initial model parameters
# obj$setParameter(logtheta  = c(init=log(theta.init),  lb=-5,    ub=ub[1]))
# obj$setParameter(mu        = c(init=mu.init,          lb=lb[2], ub=ub[2]))
# obj$setParameter(logsigmaX  = c(init=log(sigmaX.init),lb=-5,    ub=ub[3]))
# obj$setParameter(logsigmaY = c(init=log(sigmaY.init), lb=-5,    ub=ub[4]))
# 
# ctsm.opt.results = list()
# ctsm.pars.results = list()
# ctsm.cor.results = list()
# ctsm.time.results = list()
# 
# for(i in 1:N){
#   data.ctsm = data.frame(t=tsim,y=Y.kalman[i,])
#   x.init = mean(data.ctsm$y,na.rm=T) + c(0,-2,2)*sd(data.ctsm$y,na.rm=T)
#   obj$setParameter(x = c(init=x.init[1],lb=x.init[2],ub=x.init[3]))
#   out = tryCatch(
#     {
#       ctsm.time.results[[i]] = system.time(fit <- obj$estimate(data.ctsm))
#       ctsm.opt.results[[i]] = fit
#       ctsm.pars.results[[i]] = fit$xm[-1][c("logtheta","mu","logsigmaX","logsigmaY")]
#       ctsm.cor.results[[i]] = fit$corr[2:5,2:5]
#     },
#     error = function(e) print(e),
#     warning = function(w) print(w)
#   )
#   print(i)
# }

################################################################################
# COMPARE RESULTS
################################################################################

# optdata3 = list(tmb.opt.results=tmb.opt.results, tmb.pars.results=tmb.pars.results, tmb.cor.results=tmb.cor.results,
#                tmb.time.results=tmb.time.results, exact.opt.results=exact.opt.results, exact.pars.results=exact.pars.results,
#                exact.cor.results=exact.cor.results, exact.time.results=exact.time.results,ekf.opt.results=ekf.opt.results,
#                ekf.pars.results=ekf.pars.results, ekf.cor.results=ekf.cor.results, ekf.time.results=ekf.time.results,ukf.opt.results=ukf.opt.results,
#                ukf.pars.results=ukf.pars.results, ukf.cor.results=ukf.cor.results, ukf.time.results=ukf.time.results
#                )
# save(optdata3,file="optdata3.Rda")
# load("optdata3.Rda")
# for(name in names(optdata3)){
#   assign(name,optdata3[[name]])
# }

# if(!(exists("tmb.pars.results") & exists("exact.pars.results") & exists("kalman.pars.results") & exists("ctsm.pars.results"))){
#   # load("optdata.Rda")
#   load("optdata2.Rda")
#   for(name in names(optdata2)){
#     assign(name,optdata2[[name]])
#   }
# }

# true parameter values
theta  = 1
mu     = 1
sigmaX = 1e-1
sigmaY = 1e-2^2

# # boxplot for theta
# theta.data=
#   cbind(
#   tmb=sapply(tmb.pars.results,function(x) exp(x[1])),
#   exact=sapply(exact.pars.results,function(x) exp(x[1])),
#   kalman=sapply(kalman.pars.results,function(x) exp(x[1])),
#   ctsm=sapply(ctsm.pars.results,function(x) exp(x[1]))
# )
# boxplot(theta.data,main="theta")
# abline(h=theta,lw=2,col="red")
# 
# # boxplot for mu
# mu.data=
#   cbind(
#     tmb=sapply(tmb.pars.results,function(x) x[2]),
#     exact=sapply(exact.pars.results,function(x) x[2]),
#     kalman=sapply(kalman.pars.results,function(x) x[2]),
#     ctsm=sapply(ctsm.pars.results,function(x) x[2])
#   )
# boxplot(mu.data,main="mu")
# abline(h=0.5,lw=2,col="red")
# 
# # boxplot for sigmaX
# sigmaX.data=
#   cbind(
#     tmb=sapply(tmb.pars.results,function(x) exp(x[3])),
#     exact=sapply(exact.pars.results,function(x) exp(x[3])),
#     kalman=sapply(kalman.pars.results,function(x) exp(x[3])),
#     ctsm=sapply(ctsm.pars.results,function(x) exp(x[3]))
#   )
# boxplot(sigmaX.data,main="sigmaX")
# abline(h=sigmaX,lw=2,col="red")
# 
# # boxplot for sigmaY
# sigmaY.data=
#   cbind(
#     tmb=sapply(tmb.pars.results,function(x) exp(x[4])),
#     exact=sapply(exact.pars.results,function(x) exp(x[4])),
#     kalman=sapply(kalman.pars.results,function(x) exp(x[4])),
#     ctsm=sapply(ctsm.pars.results,function(x) exp(x[4]))
#   )
# boxplot(sigmaY.data,main="sigmaY")
# abline(h=sigmaY,lw=2,col="red")
# 
# time.comparison = cbind(
#   tmb = sapply(tmb.time.results,function(x) x[1]),
#   exact = sapply(exact.time.results,function(x) x[1]),
#   kalman = sapply(kalman.time.results,function(x) x[1]),
#   ctsm = sapply(ctsm.time.results,function(x) x[1])
# )
# boxplot(time.comparison,main="user time (seconds)")

#################################################################
#################### WITHOUT CTSM-R #############################
#################################################################

# parameter matrices
tmb.mat = matrix(unlist(tmb.pars.results),ncol=4,byrow=T)
exact.mat = matrix(unlist(exact.pars.results),ncol=4,byrow=T)
ekf.mat = matrix(unlist(ekf.pars.results),ncol=4,byrow=T)
ukf.mat = matrix(unlist(ukf.pars.results),ncol=4,byrow=T)
# julia.mat = matrix(unlist(julia.pars.results),ncol=4,byrow=T)

# boxplot for theta
theta.data = list(
  tmb = tmb.mat[,1],
  exact = exact.mat[,1],
  ekf = ekf.mat[,1],
  ukf = ekf.mat[,1]
  # julia = julia.mat[,1]
)
boxplot(lapply(theta.data,exp),main="theta")
abline(h=theta,lw=2,col="red")

# boxplot for mu
mu.data = list(
  tmb = tmb.mat[,2],
  exact = exact.mat[,2],
  ekf = ekf.mat[,2],
  ukf = ekf.mat[,2]
  # julia = julia.mat[,2]
)
boxplot(mu.data,main="mu",ylim=c(-1,3))
abline(h=mu,lw=2,col="red")

# boxplot for sigmaX
sigmaX.data = list(
  tmb = tmb.mat[,3],
  exact = exact.mat[,3],
  ekf = ekf.mat[,3],
  ukf = ekf.mat[,3]
  # julia = julia.mat[,3]
)
boxplot(lapply(sigmaX.data,exp),main="sigmaX")
abline(h=sigmaX,lw=2,col="red")

# boxplot for sigmaY
sigmaY.data=list(
  tmb = tmb.mat[,4],
  exact = exact.mat[,4],
  ekf = ekf.mat[,4],
  ukf = ekf.mat[,4]
  # julia = julia.mat[,4]
)
boxplot(lapply(sigmaY.data,exp),main="sigmaY")
abline(h=sigmaY,lw=2,col="red")

#boxplot for time
time.comparison = list(
 tmb = unlist(tmb.time.results),
 exact = unlist(exact.time.results),
 ekf = unlist(ekf.time.results),
 ukf = unlist(ukf.time.results)
 # julia = unlist(julia.time.results)*1000
)
boxplot(time.comparison,main="cpu_time")

rel.time.comparison = sapply(time.comparison,function(x) x/median(time.comparison$ekf))
boxplot(rel.time.comparison,main="cpu_time",ylim=c(0,100))


#boxplot for convergence status
convergence.data = cbind(
  tmb = sum(sapply(tmb.opt.results,function(x) x$message) == "false convergence (8)"),
  exact = sum(sapply(exact.opt.results,function(x) x$message)=="false convergence (8)"),
  ekf = sum(sapply(ekf.opt.results,function(x) x$message)=="false convergence (8)"),
  ukf = sum(sapply(ukf.opt.results,function(x) x$message)=="false convergence (8)")
)
barplot(convergence.data,main="convergence status")
