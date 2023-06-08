#######################################################
# FIRST FUNCTION TO RUN WHEN BUILDING
#######################################################

init_build = function(self, private) {
  
  # basic checks
  if (length(private$sys.eqs) == 0) {
    stop("Need at least one system equation")
  }
  if (length(private$obs.eqs) == 0) {
    stop("Need at least one observation equation")
  }
  if (length(private$obs.eqs)!=length(private$obs.var)) {
    stop("The number of observation equations and specified variances do not match.")
  }
  if (is.null(private$initial.state)) {
    stop("You must set an initial state estimate and covariance")
  }
  
  # set system size variables
  private$diff.processes = unique(unlist(lapply(private$sys.eqs, function(x) x$diff)))
  private$n = length(private$sys.eqs)
  private$m = length(private$obs.eqs)
  private$ng = length(private$diff.processes) - 1
  
  # add method to filepath and modelname to reflect 'method'
  private$cppfile.path.with.method = paste(private$cppfile.path,"_",private$method,sep="")
  private$modelname.with.method = paste(private$modelname,"_",private$method,sep="")
  
  return(invisible(self))
}

#######################################################
# CHECK IF ALGEBRAIC EQUATIONS ARE OK
#######################################################

check_algebraics_before_applying = function(self, private) {
  
  # So now we know all parameters, constants and inputs
  
  # An algebraic state equation
  
  # An observation can not be transformed
  
  # An observation can not depend on observations
  
  for (i in 1:length(private$alg.eqs)) {
    obj = private$alg.eqs[[i]]$form
    
    # You can not apply algebraics to a state
    if (deparse(obj[[2]]) %in% names(private$sys.eqs)) {
      stop("You can't redefine a state: ", deparse(obj))
    }
    
    # You can't apply algebraics to an input
    if (deparse(obj[[2]]) %in% names(private$inputs)) {
      stop("You can't redefine an input: ", deparse(obj))
    }
    
    # You can't apply algebraics to an observation
    if (deparse(obj[[2]]) %in% names(private$inputs)) {
      stop("You can't redefine an input: ", deparse(obj))
    }
    
  }
  return(invisible(self))
}

#######################################################
# APPLY ALGEBRAIC TRANSFORMATIONS
#######################################################

apply_algebraics = function(self, private) {
  
  # copy lists and apply algebraics
  sys.eqs = lapply(private$sys.eqs,function(x) x$rhs)
  obs.eqs = lapply(private$obs.eqs,function(x) x$rhs)
  obs.var = lapply(private$obs.var,function(x) x$rhs)
  
  algs = private$alg.eqs
  # if algs is empty just create a single fake entry to get model into
  if (is.null(algs)){
    algs = list(list(rhs=NULL))
    names(algs) = "////....i_dont_match_any_variable....////"
  }
  for (i in 1:length(algs)){
    temp_list = list()
    temp_list[[names(algs)[i]]] = algs[[i]]$rhs
    sys.eqs = lapply(sys.eqs, function(x) do.call("substitute",list(x,temp_list)))
    obs.eqs = lapply(obs.eqs, function(x) do.call("substitute",list(x,temp_list)))
    obs.var = lapply(obs.var, function(x) do.call("substitute",list(x,temp_list)))
  }
  
  # overload form and rhs, and re-add the system, observations and variances
  # systems
  sys = private$sys.eqs
  for(i in 1:length(sys.eqs)) {
    sys[[i]]$form[[3]] = sys.eqs[[i]]
    new_eq = sys[[i]]$form
    # self$add_trans_systems(new_eq)
    private$add_trans_systems(new_eq)
  }
  
  # observations
  obs = private$obs.eqs
  for(i in 1:length(obs.eqs)) {
    obs[[i]]$form[[3]] = obs.eqs[[i]]
    new_eq = obs[[i]]$form
    # self$add_trans_observations(new_eq)
    private$add_trans_observations(new_eq)
  }
  
  # observation variances
  var = private$obs.var
  for(i in 1:length(obs.var)) {
    var[[i]]$form[[3]] = obs.var[[i]]
    new_eq = var[[i]]$form
    # self$add_trans_observation_variances(new_eq)
    private$add_trans_observation_variances(new_eq)
  }
  
  return(invisible(self))
}

#######################################################
# UPDATE DIFF TERMS TO FIT ALGEBRAIC EQUATIONS
#######################################################

update_diffterms = function(self, private) {
  
  for (i in 1:length(private$sys.eqs.trans)) {
    private$diff.terms[[i]] = lapply(private$diff.processes, function(x) Deriv::Deriv(private$sys.eqs.trans[[i]]$rhs, x, cache.exp=FALSE))
    names(private$diff.terms[[i]]) = private$diff.processes
  }
  names(private$diff.terms) = names(private$sys.eqs.trans)
  
}

#######################################################
# APPLY LAMPERTI TRANSFORM IF SPECIFIED
#######################################################

apply_lamperti = function(self, private) {
  
  transform = private$lamperti$transform
  if (transform=="identity") {
    # do nothing
    return(NULL)
  }
  
  if (is.null(private$lamperti$states)){
    private$lamperti$states = private$state.names
  }
  states = private$lamperti$states
  
  # remove those states that can't have lamperti transform
  dw.bool.list = lapply(private$diff.terms, function(x) unlist(lapply(x, function(y) y!=0)))
  states_bool = rep(TRUE,length(states))
  for (i in 1:length(states)) {
    if (sum(dw.bool.list[[states[i]]]) > 2) {
      warning("Unable to perform Lamperti transformation for ",states[i], " because there are multiple diffusion terms")
      states_bool[i] = FALSE
    }
  }
  states = states[states_bool]
  dw.bool.list = dw.bool.list[states_bool]
  diff.terms = private$diff.terms[states_bool]
  
  # define and select lamperti transform from list below
  psi = list(
    log = list( quote(exp(x)), quote(1/x),quote(-1/x^2)),
    logit = list( quote(exp(x/(1+x))), quote(1/(x*(1-x))), quote((2*x-1)/(x^2*(x-1)^2))),
    'sqrt-logit' = list( quote(0.5*(sin(x)+1)), quote(1/(sqrt(x*(1-x)))), quote(0.5*(2*x-1)/(x*(1-x)*sqrt(x*(1-x)))))
    # 'power-logit' = quote(1/(x*(1-x^a)))
  )[[transform]]
  names(psi) = c("..psi..","..dpsidx..","..d2psidx2..")
  
  if (length(states)>0) {
    
    # We first need to sub in the state name in the transformation
    for (i in 1:length(states)) {
      psi = lapply(psi, function(x) do.call("substitute",list(x,list(x=parse(text=states[i])[[1]]))))
      
      dw = parse(text=names(diff.terms[[i]][dw.bool.list[[i]]])[2])[[1]]
      
      # apply lamperti
      lamperti_eq = do.call("substitute", list( quote((f * ..dpsidx.. + 0.5 * g^2 * ..d2psidx2.. ) * dt + g * ..dpsidx.. * dw),
                                                c(psi,list(f=diff.terms[[states[i]]]$dt, g=diff.terms[[states[i]]][[dw]], dw=dw))))
      # simplify causes 'dt' to be put in front of the expression
      lamperti_eq = Deriv::Simplify(lamperti_eq)
      
      # apply transformation to all instances of the state
      mylist = list(psi[['..psi..']])
      names(mylist) = states[i]
      final_eq = do.call("substitute", list(lamperti_eq,mylist))
      
      # convert to formula and parse to add_systems
      form = as.formula(paste(
        parse(text=paste("d",states[i],sep=""))[[1]],
        paste(deparse(final_eq),collapse=""),
        sep ="~"
      ))
      # self$add_trans_systems(form)
      private$add_trans_systems(form)
    }
  }
  
  return(invisible(self))
}

#######################################################
# LAST CHECK BEFORE COMPILING
#######################################################

lastcheck_before_compile = function(self, private) {
  
  # CHECK IF OBS RHS CONTAINS AT LEAST ONE STATE
  bool = unlist(lapply(private$obs.eqs.trans, function(x) any(private$state.names %in% x$allvars)))
  if (any(!bool)) {
    stop("The following observation(s) do not relate to any states: \n\t ",paste(private$obs.names[!bool],collapse=", "))
  }
  
  # CHECK IF OBS RHS DO NOT CONTAIN OBS
  bool = unlist(lapply(private$obs.eqs.trans, function(x) any(private$obs.names %in% x$allvars)))
  if (any(bool)) {
    stop("The following observation(s) attempt to observe other observations! \n\t ",paste(private$obs.names[!bool],collapse=", "))
  }
  
  # CHECK IF ALL RHS VARIABLES ARE PROVIDED BY USER
  vars = list()
  vars[[1]] = unlist(lapply(private$sys.eqs.trans, function(x) x$allvars))
  vars[[2]] = unlist(lapply(private$obs.eqs.trans, function(x) x$allvars))
  vars[[3]] = unlist(lapply(private$obs.var.trans, function(x) x$allvars))
  rhs.vars = unique(unlist(vars))
  
  # given.vars = c( private$parameter.names, private$input.names[-1], private$state.names, private$constant.names)
  given.vars = c( private$parameter.names, private$input.names, private$state.names, private$constant.names)
  bool = rhs.vars %in% given.vars
  if (any(!bool)) {
    stop("The following variables(s) are unaccounted for: \n\t ",
         paste(rhs.vars[!bool],collapse=", "))
  }
  
  # CHECK IF ANY PROVIDED VARIABLES ARE UNUSED
  given.vars = c( private$parameter.names, private$input.names[-1], private$constant.names)
  bool = given.vars %in% rhs.vars
  if (any(!bool)) {
    stop("The following variables(s) are unused: \n\t ",
         paste(given.vars[!bool],collapse=", "))
  }
  
  return(invisible(self))
}

#######################################################
# COMPILE CPP FUNCTION
#######################################################

compile_cppfile = function(self, private) {
  
  if(private$compile){
    write_cppfile(self, private)
    TMB::compile(paste(private$cppfile.path,".cpp",sep=""))
    # reload shared libraries
    try(dyn.unload(TMB::dynlib(private$cppfile.path)),silent=T)
    try(dyn.load(TMB::dynlib(private$cppfile.path)),silent=T)
  }
  
  # If compile=FALSE then the shared library must exist
  if (!private$compile) {
    
    # mac : check that files exist
    if (.Platform$OS.type=="unix") {
      modelpath.so = paste(private$cppfile.path,".so",sep="")
      if (!file.exists(modelpath.so)) {
        message("Compiling model...")
        private$compile = TRUE
        compile_cppfile(self, private)
      }
    }
    
    #windows: check that files exist
    if (.Platform$OS.type=="windows") {
      modelpath.dll = paste(private$cppfile.path,".dll",sep="")
      if (!file.exists(modelpath.dll)) {
        message("Compiling model...")
        private$compile = TRUE
        compile_cppfile(self, private)
      }
    }
    
    # reload libraries
    message("Loading pre-compiled model...")
    try(dyn.unload(TMB::dynlib(private$cppfile.path)),silent=T)
    try(dyn.load(TMB::dynlib(private$cppfile.path)),silent=T)
    
    
  }
  
  return(invisible(self))
}



