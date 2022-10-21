getModelInfo = function(model,data){
  
  sdeEq = model$sdeEq
  obsEq = model$obsEq
  obsVar = model$obsVar
  n = length(sdeEq)
  m = length(obsEq)
  
  # Get state variables
  state = c()
  for(i in 1:n){
    state[i] = str_match(deparse(sdeEq[[i]][[1]][[2]]), "d([a-zA-Z][:alnum:]*)")[2]
  }
  
  # Get right-hand sides
  rhs = list()
  for(i in 1:n){
    rhs[[i]]   = sdeEq[[i]][[1]][[3]]
  }
  
  # Get all differential processes dt, dw, dw1,...
  diff.processes = "dt"
  for(i in 1:n){
    expr = paste(deparse(rhs[[i]]),collapse="")
    diff.proc = na.omit(unlist(str_extract_all(expr, "dw([0-9]*)")))
    diff.processes = unique(c(diff.processes, diff.proc))
  }
  ng = length(diff.processes)-1
  
  # Get observation variables y,...
  obs = c()
  for(i in 1:m){
    obs[i]   = all.vars(obsEq[[i]][[1]][[2]])
  }
  
  # Get drift f and diffusion g by differentating rhs with diff.processes
  diff.terms = list()
  for(i in 1:n){
    diff.terms[[i]]        = lapply(diff.processes, function(x) { Deriv(rhs[[i]], x, cache.exp=FALSE) })
    names(diff.terms[[i]]) = diff.processes
  }
  
  returnlist = 
    list(
      sdeEq = sdeEq,
      obsEq = obsEq,
      obsVar = obsVar,
      n = n,
      m = m,
      state = state,
      rhs = rhs,
      diff.processes = diff.processes,
      ng = ng,
      obs = obs,
      diff.terms = diff.terms
      )
  
  return(returnlist)
}
