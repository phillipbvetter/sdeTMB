makeNLL = function(model,data,method="ekf",control=list(loss_function="identity",c=3),compile=TRUE,silent=TRUE,pathdir="",
                   map=list(),state.dependent.diffusion=FALSE){
  
  library(Deriv) #for Simplify command to remove e.g x/x^2 = 1/x
  library(stringr) #for easy string-matching
  
  # Substitute algebraic expressions
  for(i in 1:length(model$algeqs)){
    curlist = list()
    curlist[[names(model$algeqs)[i]]] = model$algeqs[[i]][[1]]
    model$sdeEq = lapply(model$sdeEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsEq = lapply(model$obsEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsVar = lapply(model$obsVar, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  }
  
  #Get model information (state,rhs, diff.processes, etc...)
  model_info = getModelInfo(model,data)
  names_model_info = names(model_info)
  for(i in 1:length(model_info)){
    assign(names_model_info[i],model_info[[i]])
  }
  
  # Check if model and data entries are OK
  data = checkModelAndData(model,data,control,method,state.dependent.diffusion)
  
  # MAXIMUM A POSTERIOR
  # Change provided MAPs to account for frozen parameters via tmb map-function
  if(is.null(data$MAPbool)){
    data$MAPbool = 0
  }
  if(data$MAPbool==1){
    ids = which(names(data$pars) %in% names(map))
    if(length(map)>0){
      data$MAPmean[ids] = as.vector(unlist(data$pars[ids]))
    }
  }
  
  # Loss Function
  data= applyLossFunction(data,model,control)
  
  # Overwrite modelname with path extension
  model$modelname2 = paste(pathdir,model$modelname,sep="")
  
  # Check if the system is linear in the states and observations
  IsLinear = IsSystemLinear(model,compile)
  
  ############################################################################################################
  # CASE 1 : LINEAR SYSTEM
  ############################################################################################################
  if(IsLinear){
    switch(method,
           #METHOD 1
           tmb =, #proceed to tmb_exact (same commands to be run but they are different methods!
           #METHOD 2
           tmb_exact = {
             recompFun(compile,method,model,data,control)
             tmbpars = list()
             for(i in 1:n){
               tmbpars[[state[i]]] = data[[state[i]]]
             }
             tmbpars = c(tmbpars, data$pars)
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=model$modelname, map=map, silent=silent)
             return(nll)
             # timing = system.time(opt <- nlminb(nll$par,nll$fn,nll$gr))
             # return(list(nll=nll,opt=opt,timing=timing))
           },
           #METHOD 3
           ukf = {
             recompFun(compile,method,model,data,control)
             # Prepare data
             tmbpars = data$pars
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, DLL=model$modelname, map=map, silent=silent)
             return(nll)
             # timing = system.time(opt <- nlminb(nll$par,nll$fn,nll$gr,nll$he))
             # return(list(nll=nll,opt=opt,timing=timing))
           },
           ekf = {
             recompFun(compile,method,model,data,control)
             # Prepare data
             tmbpars = data$pars
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, DLL=model$modelname, map=map, silent=silent)
             return(nll)
             # timing = system.time(opt <- nlminb(nll$par,nll$fn,nll$gr,nll$he))
             # return(list(nll=nll,opt=opt,timing=timing))
           },
           kalman_adaptive = {
             if(compile){
               write_adaptiveKalman_julia(model, data, silent=silent)
             }
             if(!compile & !julia_exists(model$modelname)){
               print("You asked me not to compile, but the julia function doesn't exist so I will compile anyway")
               write_adaptiveKalman_julia(model, data, silent=silent)
             }
             pars = unlist(data$pars)
             hyperpars = c()
             hyperpars[[1]] = c(data$X0,data$P0)
             hyperpars[[2]] = data$t
             datamat = matrix(nrow=length(data$t),ncol=m)
             for(i in 1:m){
               datamat[,i] = data[[obs[i]]]
               datamat[is.na(datamat[,1]),i] = NaN
             }
             hyperpars[[3]] = datamat
             
             timing = system.time(opt.temp <- julia_call(model$modelname,pars,hyperpars))
             opt.julia = list()
             opt.julia$par = opt.temp[[1]]
             names(opt.julia$par) = names(data$pars)
             opt.julia$objective = opt.temp[[2]]
             opt.julia$convergence = opt.temp[[3]]
             opt.julia$iterations = opt.temp[[4]]
             opt.julia$evaluations = c("function"=opt.temp[[5]],"gradient"=opt.temp[[6]],"hessian"=opt.temp[[7]])
             return(list(opt=opt.julia,timing=timing))
           }
    )
  }
  
  ############################################################################################################
  # CASE 2 : NON-LINEAR SYSTEM
  ############################################################################################################
  
  if(!IsLinear){
    if(method=="tmb_exact"){
      print("The system is non-linear so the exact method is not available. I will proceed with the approximate TMB method")
      method = "tmb"
    }
    switch(method,
           #METHOD 1
           tmb = {
             recompFun(compile,method,model,data,control)
             # Prepare data for TMB
             tmbpars = list()
             for(i in 1:n){
               tmbpars[[state[i]]] = data[[state[i]]]
             }
             tmbpars = c(tmbpars, data$pars)
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=model$modelname, map=map, silent=silent)
             return(nll)
           },
           #METHOD 2
           ukf=,
           ekf = {
             recompFun(compile,method,model,data,control)
             # Prepare data for TMB
             tmbpars = data$pars
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, DLL=model$modelname, map=map, silent=silent)
             return(nll)
           },
           kalman_adaptive = {
             if(compile){
               write_adaptiveKalman_julia(model, data, silent=silent)
             }
             pars = unlist(data$pars)
             hyperpars = c()
             hyperpars[[1]] = c(data$X0,data$P0)
             hyperpars[[2]] = data$t
             datamat = matrix(nrow=length(data$t),ncol=m)
             for(i in 1:m){
               datamat[,i] = data[[obs[i]]]
               datamat[is.na(datamat[,1]),i] = NaN
             }
             hyperpars[[3]] = datamat
             
             opt.temp = julia_call(model$modelname,pars,hyperpars)
             opt.julia = list()
             opt.julia$par = opt.temp[[1]]
             names(opt.julia$par) = names(data$pars)
             opt.julia$objective = opt.temp[[2]]
             opt.julia$convergence = opt.temp[[3]]
             opt.julia$iterations = opt.temp[[4]]
             opt.julia$evaluations = c("function"=opt.temp[[5]],"gradient"=opt.temp[[6]],"hessian"=opt.temp[[7]])
             return(opt.julia)
           }
    )
  }
}
