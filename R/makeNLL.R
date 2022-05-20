makeNLL = function(model,data,map=list(),method="TMB",exact=FALSE){
  
  # Initialize
  sdeEq = model$sdeEq
  obsEq = model$obsEq
  modelname = model$modelname
  
  state = c()
  n = length(sdeEq)
  for(i in 1:n){
    state[i] = deparse(as.name(sub("^d?([[:alnum:]]+)", "\\1", sdeEq[[i]][[1]][[2]])))
  }
  
  # C++ file name
  full_modelname = paste(modelname,".cpp",sep="")
  
  # Check if the system is linear in the states and observations or not
  IsLinear = IsSystemLinear(model)
  
  ############################################################################################################
  # CASE 1 : The system is linear and the exact matrix exponential solution is desired
  ############################################################################################################
  if(IsLinear & exact){
    write_linearExact_cpp(model,data)
    #Compile the C++ file
    compile(full_modelname)
    dyn.load(dynlib(modelname))
    # Prepare list of parameter values
    tmbpars = list()
    for(i in 1:n){
      tmbpars[[state[i]]] = data[[state[i]]]
    }
    tmbpars = c(tmbpars, data$pars)
    tmbdata = data
    # Return neg. log likelihood
    nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=modelname, map=map)
    
  }
  
  ############################################################################################################
  # CASE 2 : The system is linear but a faster solution is sought e.g. TMB-style or via Kalman filtering
  ############################################################################################################
  if(IsLinear & !exact){
    if(method=="TMB"){
      write_TMB_cpp(model,data)
      #Compile the C++ file
      compile(full_modelname)
      dyn.load(dynlib(modelname))
      # Prepare data for TMB
      tmbpars = list()
      for(i in 1:n){
        tmbpars[[state[i]]] = data[[state[i]]]
      }
      tmbpars = c(tmbpars, data$pars)
      tmbdata = data
      # Return neg. log likelihood
      nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=modelname, map=map)
    } else if (method=="kalman") {
      write_ExtendedKalman_cpp(model,data)
      #Compile the C++ file
      compile(full_modelname)
      dyn.load(dynlib(modelname))
      # Prepare data for TMB
      tmbpars = data$pars
      tmbdata = data
      # Return neg. log likelihood
      nll = MakeADFun(tmbdata, tmbpars, DLL=modelname, map=map)
    }
  }
  
  ############################################################################################################
  # CASE 3 :The system is non-linear and a solution is sought via TMB or Kalman filtering (We need other methods here)
  ############################################################################################################
  if(!IsLinear){
    if(method=="TMB"){
      write_TMB_cpp(model,data)
      #Compile the C++ file
      compile(full_modelname)
      dyn.load(dynlib(modelname))
      # Prepare data for TMB
      tmbpars = list()
      for(i in 1:n){
        tmbpars[[state[i]]] = data[[state[i]]]
      }
      tmbpars = c(tmbpars, data$pars)
      tmbdata = data
      # Return neg. log likelihood
      nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=modelname, map=map)
    } else if (method=="kalman"){
      write_ExtendedKalman_cpp(model,data)
      #Compile the C++ file
      compile(full_modelname)
      dyn.load(dynlib(modelname))
      # Prepare data for TMB
      tmbpars = data$pars
      tmbdata = data
      # Return neg. log likelihood
      nll = MakeADFun(tmbdata, tmbpars, DLL=modelname, map=map)
      }
  }
  
  
  return(nll)
}
