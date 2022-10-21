# Initialize
rm(list=ls())

# Libraries
library(devtools)
library(Deriv)
library(TMB)
library(bench)

# Set directory
this_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(this_path))

# Loading functions package style 
pack = as.package("../R")
load_all(pack)

# Construct model list
model = list()

model$sdeEq = list()
model$sdeEq[[1]] = expression( dx1 ~ x2 * dt + sigma_x1 * dw1 )
model$sdeEq[[2]] = expression( dx2 ~ (lambda * (1-x1^2)*x2 - x1) * dt + sigma_x2 * dw2 )
# Note: Currently not allowed to use dw, must be numbered e.g. dw1, dw2,...

model$obsEq = list()
model$obsEq[[1]] = expression(y1 ~ x1)
model$obsEq[[2]] = expression(y2 ~ x2)

model$obsVar = list()
model$obsVar[[1]] = expression(sigma_y1)
model$obsVar[[2]] = expression(sigma_y2)

model$algeqs = list()
model$algeqs[["lambda"]] = expression(exp(loglambda))
model$algeqs[["sigma_x1"]] = expression(exp(logsigma_x1))
model$algeqs[["sigma_x2"]] = expression(exp(logsigma_x2))
model$algeqs[["sigma_y2"]] = expression(exp(logsigma_y1))
model$algeqs[["sigma_y1"]] = expression(exp(logsigma_y2))

rBM <- function(t) rnorm(length(t),mean=0,sd=sqrt(diff(c(0,t))))

euler <- function(f=function(x) 0 , g=function(x) 1 ,tvec, X0)
{
  ## Simulation duration and time step
  Nt <- length(tvec)
  Np <- length(X0)
  
  ## Initialize horiz. and vert. position
  X <- array(NA,c(Nt,Np))
  X[1,] <- X0
  
  dB = X
  for(i in 1:Np) dB[,i] = rBM(tvec)
  
  
  ## Main time loop
  for(i in 2:Nt)
  {
    ## Update vertical position; reflect at bottom
    X[i,] <- X[i-1,] + f(X[i-1,]) * (tvec[i]-tvec[i-1]) + g(X[i-1,]) * dB[i,]
  }
  return(X)
}

# time variables
Tstep  = 1e-2   
Tfinal = 20      
Tobs   = 0.5

# Time vector
tsim <- seq(0, Tfinal, Tstep)

# Measurements with sample interval Tobs
iobs <- seq(1, length(tsim), round(Tobs/Tstep))

# parameters
lambda = 2
sigma_x1 = 0.1
sigma_x2 = 0.5
sigma_obs = 0.1

f = function(x) c( x[2], lambda * (1-x[1]^2) * x[2] - x[1] )
g = function(x) c( sigma_x1, sigma_x2 )
h = function(x) c(x[1],x[2])


# Simulate trajectories
X0 = c(1,1)
N = 10
X = matrix(nrow=N*2,ncol=length(tsim))
Y.data = matrix(NA,nrow=N*2,ncol=length(tsim))
for(i in 1:N) {
  xsim = euler(f,g,tsim,X0)
  X[(2*i-1),] = xsim[,1]
  X[(2*i),] = xsim[,2]
  Y.data[(2*i-1),iobs] <- rnorm(length(iobs), mean=X[(2*i-1),iobs], sd = sigma_obs)
  Y.data[(2*i),iobs] <- rnorm(length(iobs), mean=X[2*i,iobs], sd = sigma_obs)
}

plot(tsim,X[1,],type="l",lwd=2,lty="solid",ylim=c(-5,5))
lines(tsim,X[2,],type="l",lwd=2,lty="solid",ylim=c(-5,5))
points(tsim,Y.data[1,],col="blue",cex=1,pch=16)
points(tsim,Y.data[2,],col="red",cex=1,pch=16)

################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################
################ TMB & EXACT DATA ################

data.tmb = list()
data.tmb$t = tsim
data.tmb$y1 = Y.data[1,]
data.tmb$y2 = Y.data[2,]
data.tmb$x1 = rep(mean(data.tmb$y1,na.rm=T),length(data.tmb$y1))
data.tmb$x2 = rep(mean(data.tmb$y2,na.rm=T),length(data.tmb$y2))
data.tmb$pars = list(
  loglambda = log(2),
  logsigma_x1 = log(0.25),
  logsigma_x2 = log(0.25),
  logsigma_y1 = log(0.2),
  logsigma_y2 = log(0.2)
)

################ KALMAN DATA ################
################ KALMAN DATA ################
################ KALMAN DATA ################

data.kalman = list()
data.kalman$t = tsim
data.kalman$dt = diff(tsim)[1]/5
data.kalman$y1 = Y.data[1,]
data.kalman$y2 = Y.data[2,]
data.kalman$X0 = c( mean(data.kalman$y1,na.rm=T) , mean(data.kalman$y2,na.rm=T) )
data.kalman$P0 = cov(t(Y.data[1:2,]),use="pairwise.complete.obs")
data.kalman$pars = list(
  loglambda = log(2),
  logsigma_x1 = log(0.25),
  logsigma_x2 = log(0.25),
  logsigma_y1 = log(0.2),
  logsigma_y2 = log(0.2)
)

# optimization upper and lower bounds
lb = c( log(1e-1), log(1e-10), log(1e-10), log(1e-10), log(1e-10) )
ub = c( log(10), log(2), log(2), log(2), log(2) )

################################################################################
# CASE 1 - RESULTS USING TMB-STYLE
################################################################################

model$modelname = "nll_tmb_vdp"
# nll = makeNLL(model,data.tmb,method="tmb",compile=TRUE,silent=TRUE)

tmb.opt.results = list()
tmb.pars.results = list()
tmb.cor.results = list()
tmb.time.results = list()

for(i in 1:N){
  data.tmb$y1 = Y.data[(2*i-1),]
  data.tmb$y2 = Y.data[2*i,]
  data.tmb$x1 = rep(mean(data.tmb$y1,na.rm=T),length(data.tmb$y1))
  data.tmb$x2 = rep(mean(data.tmb$y2,na.rm=T),length(data.tmb$y2))
  # Return neg. log likelihood
  nll = makeNLL(model,data.tmb,method="tmb",compile=FALSE,silent=TRUE)
  out = tryCatch(
    {
      tmb.time.results[[i]] = mark(opt <- nlminb(nll$par,nll$fn,nll$gr,lower=lb,upper=ub),iterations=1)$median
      tmb.opt.results[[i]] = opt
      tmb.pars.results[[i]] = opt$par
      sdr = sdreport(nll)
      tmb.cor.results[[i]] = cov2cor(sdr$cov.fixed)
    },
    error = function(e) print(e)
  )
  print(i)
}

################################################################################
# CASE 3 - RESULTS USING KALMAN FILTER-STYLE
################################################################################

model$modelname = "nll_ekf_vdp"
# nll = makeNLL(model,data.kalman,method="ekf",compile=TRUE,silent=TRUE)

ekf.pars.results = list()
ekf.cor.results = list()
ekf.opt.results = list()
ekf.time.results = list()

for(i in 1:N){
  data.kalman$y1 = Y.data[(2*i-1),]
  data.kalman$y2 = Y.data[2*i,]
  data.kalman$X0 = c( mean(data.kalman$y1,na.rm=T), mean(data.kalman$y2,na.rm=T) )
  data.kalman$P0 = cov(cbind(data.kalman$y1,data.kalman$y2),use="pairwise.complete.obs")
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
    error = function(e) print(e)
    # warning = function(w) {print(w);print(solve(nll$he(opt$par)))}
  )
  print(i)
}

################################################################################
# CASE 3 - RESULTS USING UNSCENTED KALMAN FILTER-STYLE
################################################################################

model$modelname = "nll_ukf_vdp"
# nll = makeNLL(model,data.kalman,method="ukf",compile=TRUE,silent=TRUE)

ukf.pars.results = list()
ukf.cor.results = list()
ukf.opt.results = list()
ukf.time.results = list()

for(i in 1:N){
  data.kalman$y1 = Y.data[(2*i-1),]
  data.kalman$y2 = Y.data[2*i,]
  data.kalman$X0 = c( mean(data.kalman$y1,na.rm=T), mean(data.kalman$y2,na.rm=T) )
  data.kalman$P0 = cov(cbind(data.kalman$y1,data.kalman$y2),use="pairwise.complete.obs")
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
    error = function(e) print(e)
    # warning = function(w) {print(w);print(solve(nll$he(opt$par)))}
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

# parameter matrices
tmb.mat = matrix(unlist(tmb.pars.results),ncol=5,byrow=T)
ekf.mat = matrix(unlist(ekf.pars.results),ncol=5,byrow=T)
ukf.mat = matrix(unlist(ukf.pars.results),ncol=5,byrow=T)
# julia.mat = matrix(unlist(julia.pars.results),ncol=4,byrow=T)

# boxplot for lambda
lambda.data = list(
  tmb = tmb.mat[,1],
  ekf = ekf.mat[,1],
  ukf = ekf.mat[,1]
  # julia = julia.mat[,1]
)
boxplot(lapply(lambda.data,exp),main="lambda")
abline(h=lambda,lw=2,col="red")

# boxplot for sigmax1
sigmax1.data = list(
  tmb = tmb.mat[,2],
  ekf = ekf.mat[,2],
  ukf = ukf.mat[,2]
  # julia = julia.mat[,3]
)
boxplot(lapply(sigmax1.data,exp),main="sigmax1")
abline(h=sigma_x1,lw=2,col="red")

# boxplot for sigmax1
sigmax2.data = list(
  tmb = tmb.mat[,3],
  ekf = ekf.mat[,3],
  ukf = ukf.mat[,3]
  # julia = julia.mat[,3]
)
boxplot(lapply(sigmax2.data,exp),main="sigmax2")
abline(h=sigma_x2,lw=2,col="red")

# boxplot for sigmay1
sigmay1.data=list(
  tmb = tmb.mat[,4],
  ekf = ekf.mat[,4],
  ukf = ukf.mat[,4]
  # julia = julia.mat[,4]
)
boxplot(lapply(sigmay1.data,exp),main="sigma_y1")
abline(h=sigma_obs,lw=2,col="red")

# boxplot for sigmay2
sigmay2.data=list(
  tmb = tmb.mat[,5],
  ekf = ekf.mat[,5],
  ukf = ukf.mat[,5]
  # julia = julia.mat[,4]
)
boxplot(lapply(sigmay2.data,exp),main="sigma_y2")
abline(h=sigma_obs,lw=2,col="red")

#boxplot for time
time.comparison = list(
  tmb = unlist(tmb.time.results),
  ekf = unlist(ekf.time.results),
  ukf = unlist(ukf.time.results)
)
boxplot(time.comparison,main="cpu_time",ylim=c(0,4))

rel.time.comparison = sapply(time.comparison,function(x) x/median(time.comparison$ekf))
boxplot(rel.time.comparison,main="cpu_time",ylim=c(0,100))


#boxplot for convergence status
convergence.data = cbind(
  tmb = sum(sapply(tmb.opt.results,function(x) x$message) == "false convergence (8)"),
  ekf = sum(sapply(ekf.opt.results,function(x) x$message)=="false convergence (8)"),
  ukf = sum(sapply(ukf.opt.results,function(x) x$message)=="false convergence (8)")
)
barplot(convergence.data,main="convergence status")

