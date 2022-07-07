IsSystemLinear = function(model,silent=FALSE){
  
  sdeEq = model$sdeEq
  obsEq = model$obsEq
  
  n = length(sdeEq)
  m = length(obsEq)
  # Get state and observation variables
  states = c()
  rhs = list()
  for(i in 1:n){
    states[i] = deparse(as.name(sub("^d?([[:alnum:]]+)", "\\1", sdeEq[[i]][[1]][[2]])))
    rhs[[i]]   = sdeEq[[i]][[1]][[3]]
  }
  # Get the drift and diffusion terms (dt, dw1, dw2...)
  diff.processes = c("dt",sprintf(rep("dw%i",n),1:n))
  diff.terms = list()
  for(i in 1:n){
    diff.terms[[i]]        = lapply(diff.processes, function(x) { D(rhs[[i]], x) })
    names(diff.terms[[i]]) = diff.processes
  }
  # Get the observation right hand side
  obs.terms = lapply(obsEq,function(x) x[[1]][[3]])
  
  # We first check the stochastic differential equations
  for(i in 1:n){
    # 1. Drift terms
    # Select Eq
    eq       = diff.terms[[i]]$dt
    prev.vars = all.vars(eq)
    present.states = states[states %in% prev.vars]
    # Diff wrt all states
    A        = lapply(states,function(x) D(eq,x) )
    names(A) = states
    B = lapply(A,all.vars)
    # Cheeck if previously present state were lost after differentiation
    C = all(unlist(lapply(present.states,function(x) !any(x==B[[x]]))))
    if(!C){
      if(!silent){print("The system is nonlinear")}
      return(linear=FALSE)
    }
    # 2. Diffusion terms
    for(j in 1:n){
      eq = diff.terms[[i]][[j+1]]
      prev.vars = all.vars(eq)
      present.states = states[states %in% prev.vars]
      if( length(present.states) > 0 ){
        if(!silent){print("The system is nonlinear")}
        return(linear=FALSE)
      }
    }
  }
  # We secondly check the observation equations
  for(i in 1:m){
    eq = obs.terms[[i]]
    prev.vars = all.vars(eq)
    present.states = states[states %in% prev.vars]
    A        = lapply(states,function(x) D(eq,x) )
    names(A) = states
    B = lapply(A,all.vars)
    # Cheeck if previously present states were lost after differentiation
    C = all(unlist(lapply(present.states,function(x) !any(x==B[[x]]))))
    if(!C){
      if(!silent){print("The system is nonlinear")}
      return(linear=FALSE)
    }
  }
  
  if(!silent){print("The system is linear")}
  return(linear=TRUE)
}
