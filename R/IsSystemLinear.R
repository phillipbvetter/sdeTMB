IsSystemLinear = function(model,compile){
  
  returnString = c()
  linearFlag = TRUE
  
  silent = FALSE
  if(!compile){
    silent = TRUE
  }
  
  #Get info
  model_info = getModelInfo(model,data)
  names_model_info = names(model_info)
  for(i in 1:length(model_info)){
    assign(names_model_info[i],model_info[[i]])
  }
  
  # We first check the stochastic differential equations
  for(i in 1:n){
    # 1. Drift terms
    eq       = diff.terms[[i]]$dt
    prev.vars = all.vars(eq)
    present.states = state[state %in% prev.vars]
    # Derivative with respect to states
    A        = lapply(state,function(x) D(eq,x) )
    names(A) = state
    B = lapply(A,all.vars)
    # If the state is not lost after differentiation, it is a nonlinear system
    C = all(unlist(lapply(present.states,function(x) !any(x==B[[x]]))))
    if(!C){
      returnString = append(returnString, "The drift function is nonlinear")
      linearFlag = FALSE
    }
  }
  # 2. Diffusion terms
  # An SDE is considered non-linear if the diffusion is state dependent (no derivative needs to be taken)
  n.diff = length(diff.processes)-1
  for(i in 1:n){
    for(j in 1:n.diff){
      eq = diff.terms[[i]][[j+1]]
      prev.vars = all.vars(eq)
      present.states = state[state %in% prev.vars]
      if( length(present.states) > 0 ){
        returnString = append(returnString, "The diffusion function is nonlinear")
        linearFlag = FALSE
      }
    }
  } 
  
  # Observation Equations
  for(i in 1:m){
    eq = obsEq[[i]][[1]][[3]]
    prev.vars = all.vars(eq)
    present.states = state[state %in% prev.vars]
    A        = lapply(state,function(x) D(eq,x) )
    names(A) = state
    B = lapply(A,all.vars)
    # Cheeck if previously present states were lost after differentiation
    C = all(unlist(lapply(present.states,function(x) !any(x==B[[x]]))))
    if(!C){
      returnString = append(returnString, "The observation equation is nonlinear")
      linearFlag = FALSE
    }
  }
  
  if(linearFlag){
    if(!silent){
      message("SYSTEM IS LINEAR")
    }
    return(linear=TRUE)
  } else {
    if(!silent){
      message("SYSTEM IS NONLINEAR:")
      message(paste(returnString,collapse="\n"))
    }
    return(linear=FALSE)
  }
}
