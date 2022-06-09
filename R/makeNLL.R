makeNLL = function(model,data,map=list(),method="TMB",compile=TRUE,silent=TRUE){
  # method are "TMB", "kalman" and "exactTMB"
  
  # Initialize
  sdeEq = model$sdeEq
  obsEq = model$obsEq
  modelname = model$modelname

  state = c()
  n = length(sdeEq)
  for(i in 1:n){
    state[i] = deparse(as.name(sub("^d?([[:alnum:]]+)", "\\1", sdeEq[[i]][[1]][[2]])))
  }
  
  # Check if the system is linear in the states and observations or not
  IsLinear = IsSystemLinear(model)
  
  # Function for recompilation
  recompFun = function(compile,method,model,data,modelname){
    full_modelname = paste(modelname,".cpp",sep="")
    if(compile){
      switch(method,
             "TMB" = write_TMB_cpp(model,data),
             "kalman" = write_ExtendedKalman_cpp(model,data),
             "TMBexact" = write_linearExact_cpp(model,data)
      )
      compile(full_modelname)
      dyn.load(dynlib(modelname))
    }
    if(!file.exists(full_modelname)){
      print("You asked me not to compile, but the model C++ file doesn't exist so I will compile anyway")
      recompFun(TRUE,method,model,data,modelname)
    }
    if(!any(modelname==names(getLoadedDLLs()))) dyn.load(dynlib(modelname))
  }
  
  ############################################################################################################
  # CASE 1 : The system is linear and the exact matrix exponential solution is desired
  ############################################################################################################
if(IsLinear){
  if(method=="TMBexact"){
    recompFun(compile,method,model,data,modelname)
    # dyn.load(dynlib(modelname))
    # Prepare list of parameter values
    tmbpars = list()
    for(i in 1:n){
      tmbpars[[state[i]]] = data[[state[i]]]
    }
    tmbpars = c(tmbpars, data$pars)
    tmbdata = data
    # Return neg. log likelihood
    nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=modelname, map=map,silent=silent)
  }
  if(method=="TMB"){
    recompFun(compile,method,model,data,modelname)
    # dyn.load(dynlib(modelname))
    # Prepare data for TMB
    tmbpars = list()
    for(i in 1:n){
      tmbpars[[state[i]]] = data[[state[i]]]
    }
    tmbpars = c(tmbpars, data$pars)
    tmbdata = data
    # Return neg. log likelihood
    nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=modelname, map=map,silent=silent)
  }
  if (method=="kalman") {
    recompFun(compile,method,model,data,modelname)
    # dyn.load(dynlib(modelname))
    # Prepare data for TMB
    tmbpars = data$pars
    tmbdata = data
    # Return neg. log likelihood
    nll = MakeADFun(tmbdata, tmbpars, DLL=modelname, map=map,silent=silent)
  }
}

  ############################################################################################################
  # CASE 3 :The system is non-linear and a solution is sought via TMB or Kalman filtering (We need other methods here)
  ############################################################################################################
  if(!IsLinear){
    if(method=="TMBexact"){
      print("The system is non-linear so the exact method is not available. I will proceed with the approximate TMB method")
      method = "TMB"
    }
    if(method=="TMB"){
      recompFun(compile,method,model,data,modelname)
      # dyn.load(dynlib(modelname))
      # Prepare data for TMB
      tmbpars = list()
      for(i in 1:n){
        tmbpars[[state[i]]] = data[[state[i]]]
      }
      tmbpars = c(tmbpars, data$pars)
      tmbdata = data
      # Return neg. log likelihood
      nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=modelname, map=map,silent=silent)
    }
    if(method=="kalman"){
      recompFun(compile,method,model,data,modelname)
      # dyn.load(dynlib(modelname))
      # Prepare data for TMB
      tmbpars = data$pars
      tmbdata = data
      # Return neg. log likelihood
      nll = MakeADFun(tmbdata, tmbpars, DLL=modelname, map=map,silent=silent)
      }
  }
  
  return(nll)
}
