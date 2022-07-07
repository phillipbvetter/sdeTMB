makeNLL = function(model,data,map=list(),method="TMB",compile=TRUE,silent=TRUE){
  # method are "TMB", "kalman" and "exactTMB"
  
  # Substitute algebraic expressions
  for(i in 1:length(model$algeqs)){
    curlist = list()
    curlist[[names(model$algeqs)[i]]] = model$algeqs[[i]][[1]]
    model$sdeEq = lapply(model$sdeEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsEq = lapply(model$obsEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsVar = lapply(model$obsVar, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  }
  
  # Initialize
  sdeEq = model$sdeEq
  obsEq = model$obsEq

  state = c()
  n = length(sdeEq)
  for(i in 1:n){
    state[i] = deparse(as.name(sub("^d?([[:alnum:]]+)", "\\1", sdeEq[[i]][[1]][[2]])))
  }
  
  # Check if the system is linear in the states and observations or not
  Shhh = FALSE
  if(!compile){
    Shhh = TRUE
  }
  IsLinear = IsSystemLinear(model,Shhh)
  
  # Function for recompilation
  recompFun = function(compile,method,model,data){
    full_modelname = paste(model$modelname,".cpp",sep="")
    if(compile){
      switch(method,
             TMB = write_TMB_cpp(model,data),
             kalman = write_ExtendedKalman_cpp(model,data),
             TMBexact = write_linearExact_cpp(model,data)
      )
      compile(full_modelname)
      #reload the library
      try(dyn.unload(dynlib(model$modelname)),silent=T)
      dyn.load(dynlib(model$modelname))
    }
    if(!file.exists(full_modelname)){
      print("You asked me not to compile, but the model C++ file doesn't exist so I will compile anyway")
      recompFun(TRUE,method,model,data)
    }
    #if you restart R and do not require recompilation, then just load the library
    if(!any(model$modelname==names(getLoadedDLLs()))) dyn.load(dynlib(model$modelname))
  }
  
  ############################################################################################################
  # CASE 1 : The system is linear and the exact matrix exponential solution is desired
  ############################################################################################################
if(IsLinear){
  switch(method,
         #METHOD 1
         TMB =, #TMB is identical to TMBexact
         #METHOD 2
         TMBexact = {
           recompFun(compile,method,model,data)
           tmbpars = list()
           for(i in 1:n){
             tmbpars[[state[i]]] = data[[state[i]]]
           }
           tmbpars = c(tmbpars, data$pars)
           tmbdata = c(data,data$constants)
           # Return neg. log likelihood
           nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=model$modelname, map=map, silent=silent)
         },
         #METHOD 3
         kalman = {
           recompFun(compile,method,model,data)
           # Prepare data
           tmbpars = data$pars
           tmbdata = c(data,data$constants)
           # Return neg. log likelihood
           nll = MakeADFun(tmbdata, tmbpars, DLL=model$modelname, map=map, silent=silent)
         }
  )
}
  
  ############################################################################################################
  # CASE 2 :The system is non-linear and a solution is sought via TMB or Kalman filtering (We need other methods here)
  ############################################################################################################
  
  if(!IsLinear){
    if(method=="TMBexact"){
      print("The system is non-linear so the exact method is not available. I will proceed with the approximate TMB method")
      method = "TMB"
    }
    switch(method,
           #METHOD 1
           TMB = {
             recompFun(compile,method,model,data)
             # Prepare data for TMB
             tmbpars = list()
             for(i in 1:n){
               tmbpars[[state[i]]] = data[[state[i]]]
             }
             tmbpars = c(tmbpars, data$pars)
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, random=state, DLL=model$modelname, map=map, silent=silent)
           },
           #METHOD 2
           kalman = {
             recompFun(compile,method,model,data,modelname)
             # Prepare data for TMB
             tmbpars = data$pars
             tmbdata = c(data,data$constants)
             # Return neg. log likelihood
             nll = MakeADFun(tmbdata, tmbpars, DLL=modelname, map=map, silent=silent)
           }
    )
  }
  return(nll)
}
