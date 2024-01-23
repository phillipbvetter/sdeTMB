# These data checking and setting functions are called automatically when the user
# requests an estimation, prediction or construct_nll.

#######################################################
# CHECK AND SET DATA BEFORE OPTIMIZATION
#######################################################

basic_data_check = function(data, self, private) {
  
  # check if data is list or data.frame
  if (!(is.list(data) | is.data.frame(data))) {
    stop("The data should be a data.frame or a list")
  }
  
  # check if all inputs are in the data
  bool = private$input.names %in% names(data)
  if (any(!bool)){
    stop("The following required inputs were not provided in the data: 
         ", paste(private$input.names[!bool],collapse=", "))
  }
  
  # check if all (basic) observations are in the data.
  # The basic observations are not obs.names, but variables on lhs of obs equations, 
  # due to the ability to say e.g. log(y) ~ ..., obsname=log(y)
  required.obs = unique(unlist(sapply(private$obs.eqs.trans, function(ls) all.vars(ls$lhs))))
  bool = required.obs %in% names(data)
  if (any(!bool)){
    stop("The following required observations were not provided in the data: 
         ", required.obs[!bool])
  }
  
  # time vector must be increasing
  if (any(diff(data$t)<=0)) {
    ids = which(diff(data$t)<=0)
    stop(sprintf("The time-vector is non-increasing at the following indice(s) %s",paste(ids,collapse=", ")))
  }
  
  # time vector must only contain numerics
  if (any(is.na(data$t))) {
    ids = which(is.na(data$t))
    stop(sprintf("The time-vector is NA at the following indice(s) %s",paste(ids,collapse=", ")))
  }
  
  return(invisible(self))
}

#######################################################
# CALCULATE THE COMPLEX OBSERVATION NAMES
#######################################################

calculate_complex_observation_lefthandsides = function(data, self, private){
  
  # The class of the quote(log(x)) is 'call' whereas quote(x) is 'name', so complex
  # observation equations can be identified as being class 'call'
  
  # detect the complex obs lhs
  bool = as.vector(unlist(lapply(private$obs.eqs.trans, function(ls) inherits(ls$lhs, "call"))))
  
  # if there are none, return
  if(!any(bool)){
    return(data)
  }
  
  # otherwise, calculate these variables using data variables
  temp.data = list()
  for(i in seq_along(private$obs.eqs.trans)[bool]){
    
    # get name and lhs
    lhs = private$obs.eqs.trans[[i]]$lhs
    name = private$obs.eqs.trans[[i]]$name
    
    # Check if the variables are available in data
    bool = all.vars(lhs) %in% names(data)
    if(!all(bool)){
      stop("Unable to compute the observation ",name," because the following variable(s) are not in the provided data: 
           ",all.vars(lhs)[!bool])
    }
    
    # Compute the complex observation
    new.data.entry = with(data, eval(lhs))
    
    # append to the data
    temp.data[[name]] = new.data.entry
  }
  
  # concatenate data frames
  newdata = data.frame(data, temp.data)
  
  # return
  return(newdata)
}

#######################################################
# SETTINGS FOR ODE TIMESTEP 
#######################################################

set_ode_timestep = function(data, self, private){
  
  # If the required number of steps is N + epsilon or larger (e.g. 3+0.01) then increase step by 1, and reduce timestep there.
  # :::::EXAMPLE:::::
  # data$t = [0 , 1 , 2, 4.5], so data.dt = [1,1,2.5]
  # ode.timestep = 1. There are therefore [1, 1, 2.5] steps required. The last (2.5) has residual larger than 0.05 so (2.5 %% 1 = 0.5 > 0.05)
  # so we round up the number of steps there i.e. ode.N = [1 , 1 , 3]. The last entry is the important one. 
  # We take 3 steps, so for last entry, we must reduce the step-size to data.dt[3] / ode.N[3] = 2.5 / 3  = 0.88883333 
  # down from the set ode.timestep = 1
  
  
  # check for no provided time.step
  if(is.null(private$ode.timestep)){
    
    private$ode.timestep = diff(data$t)
    
  }
  
  # check that ode.timestep has length 1 or at least nrow(data)-1.
  if (length(private$ode.timestep) == 1) {
    
    # Recycle to correct length
    private$ode.timestep = rep(private$ode.timestep, nrow(data)-1)
    
  }  else if (length(private$ode.timestep) == nrow(data) -1) {
    
    # do nothing
    
  } else if (length(private$ode.timestep) > nrow(data) - 1 ) {
    
    # trim if it is larger than nrow(data)-1
    private$ode.timestep = head(private$ode.timestep, nrow(data)-1)
    warning("The provided ode.timestep was longer than nrow(data) - 1, only using first nrow(data)-1 entries.")
    
  } else {
    
    # throw error
    stop("Error: The provided ode.timestep must have length 1 or nrow(data)-1")
  }
  
  
  # Data time-difference, and ode.timesteps (number of steps to take)
  data.dt = diff(data$t)
  ode.timesteps = rep(1,length(data.dt))
  
  # For the data time-gaps larger than ode.timestep, we set the time-step to the requested ode.timestep
  # For all others just take a timestep equal to the time-gap in the data (which will be smaller than ode.timestep)
  ode_timestep_size = data.dt
  bool = data.dt > private$ode.timestep
  ode_timestep_size[bool] = private$ode.timestep[bool]
  
  # where data.dt > ode.timestep, calculate the required number of steps using ode.timestep
  ode.timesteps[bool] = data.dt[bool] / private$ode.timestep[bool]
  
  # Find those indices where the number of steps must increase
  epsilon.step  = 1e-2
  residual.step.bool = (ode.timesteps %% 1) > epsilon.step 
  
  # increment these steps from N to N+1
  ode.timesteps[residual.step.bool] = ceiling(ode.timesteps[residual.step.bool]) 
  
  # round other steps down (e.g. any number less than N+epsilon becomes N)
  ode.timesteps[!residual.step.bool] = floor(ode.timesteps[!residual.step.bool])
  
  # now reduce the timestep at these locations to match the integer number of steps such that ode.timesteps * ode.timestep = data.dt
  # e.g. ode.timestep = data.dt / ode.timesteps
  ode_timestep_size[residual.step.bool] = data.dt[residual.step.bool] / ode.timesteps[residual.step.bool]
  
  # store the calculated step-size, number of steps, and cumulative number of steps
  private$ode.timestep.size = ode_timestep_size
  private$ode.timesteps = ode.timesteps
  private$ode.timesteps.cumsum = c(0,cumsum(private$ode.timesteps)) #this is used in the laplace method
  
  # return
  return(invisible(self))
}


#######################################################
# IOBS VECTOR FOR LAPLACE TO AVOID USING IS.NA
#######################################################

# The reason we want to avoid using is.na is probably that it enables us
# to use the one-step-residual function from TMB...???

set_data_for_laplace_method = function(data, self, private){
  
  # Create vector to avoid using 'is.na' in TMB objective function in C++
  iobs = list()
  for (i in seq_along(private$obs.names)) {
    iobs[[i]] = seq_along(data$t)[!is.na(data[[private$obs.names[i]]])] - 1 # -1 because C++ 0-indexed
  }
  names(iobs) = paste("iobs_",private$obs.names,sep="")
  private$iobs = iobs
  
  # Did the user supply random effects initial values in the data?
  # Otherwise construct from the set_initial value
  bool = !(private$state.names %in% names(data))
  if(private$method=="laplace"){
    
    if(any(bool)){
      data[private$state.names[bool]] =
        matrix(rep(private$initial.state$x0[bool], times=length(data$t)),ncol=sum(bool))
    }
    
    # save state values as tmb initial state values
    for(i in seq_along(private$state.names)){
      
      private$tmb.initial.state.for.parameters[[i]] =
        rep(data[[private$state.names[i]]], times = c(private$ode.timesteps,1))
      
    }
    names(private$tmb.initial.state.for.parameters) = private$state.names
    
  }
  
  # return
  return(invisible(self))
}

#######################################################
# SETTINGS FOR ODE TIMESTEP 
#######################################################

set_simulation_timestep = function(data, self, private){
  
  # If the required number of steps is N + epsilon or larger (e.g. 3+0.01) then increase step by 1, and reduce timestep there.
  # :::::EXAMPLE:::::
  # data$t = [0 , 1 , 2, 4.5], so data.dt = [1,1,2.5]
  # ode.timestep = 1. There are therefore [1, 1, 2.5] steps required. The last (2.5) has residual larger than 0.05 so (2.5 %% 1 = 0.5 > 0.05)
  # so we round up the number of steps there i.e. ode.N = [1 , 1 , 3]. The last entry is the important one. 
  # We take 3 steps, so for last entry, we must reduce the step-size to data.dt[3] / ode.N[3] = 2.5 / 3  = 0.88883333 
  # down from the set ode.timestep = 1
  
  
  # check for no provided time.step
  if(is.null(private$simulation.timestep)){
    private$simulation.timestep = diff(data$t)
    
  }
  
  # check that simulation.timestep has length 1 or at least nrow(data)-1.
  if (length(private$simulation.timestep) == 1) {
    
    # Recycle to correct length
    private$simulation.timestep = rep(private$simulation.timestep, nrow(data)-1)
    
  }  else if (length(private$simulation.timestep) == nrow(data) -1) {
    
    # do nothing
    
  } else if (length(private$simulation.timestep) > nrow(data) - 1 ) {
    
    # trim if it is larger than nrow(data)-1
    private$simulation.timestep = head(private$simulation.timestep, nrow(data)-1)
    warning("The provided simulation.timestep was longer than nrow(data) - 1, only using first nrow(data)-1 entries.")
    
  } else {
    
    # throw error
    stop("Error: The provided simulation.timestep must have length 1 or nrow(data)-1")
  }
  
  
  # Data time-difference, and simulation.timesteps (number of steps to take)
  data.dt = diff(data$t)
  simulation.timesteps = rep(1,length(data.dt))
  
  # For the data time-gaps larger than ode.timestep, we set the time-step to the requested ode.timestep
  # For all others just take a timestep equal to the time-gap in the data (which will be smaller than ode.timestep)
  simulation_timestep_size = data.dt
  bool = data.dt > private$simulation.timestep
  simulation_timestep_size[bool] = private$simulation.timestep[bool]
  
  # where data.dt > ode.timestep, calculate the required number of steps using ode.timestep
  simulation.timesteps[bool] = data.dt[bool] / private$simulation.timestep[bool]
  
  # Find those indices where the number of steps must increase
  epsilon.step  = 1e-2
  residual.step.bool = (simulation.timesteps %% 1) > epsilon.step 
  
  # increment these steps from N to N+1
  simulation.timesteps[residual.step.bool] = ceiling(simulation.timesteps[residual.step.bool]) 
  
  # round other steps down (e.g. any number less than N+epsilon becomes N)
  simulation.timesteps[!residual.step.bool] = floor(simulation.timesteps[!residual.step.bool])
  
  # now reduce the timestep at these locations to match the integer number of steps such that simulation.timesteps * ode.timestep = data.dt
  # e.g. ode.timestep = data.dt / simulation.timesteps
  simulation_timestep_size[residual.step.bool] = data.dt[residual.step.bool] / simulation.timesteps[residual.step.bool]
  
  # store the calculated step-size, number of steps, and cumulative number of steps
  private$simulation.timestep.size = simulation_timestep_size
  private$simulation.timesteps = simulation.timesteps
  
  # return
  return(invisible(self))
}

set_ukf_parameters = function(unscented_hyperpars, self, private)
{
  if(is.null(unscented_hyperpars))
  {
    unscented_hyperpars = list(alpha=1, beta=0, kappa=3-private$number.of.states)
  }
  private$set_ukf_hyperpars(unscented_hyperpars)
  
  return(invisible(self))
}


set_pred_initial_state = function(initial.state, self, private)
{
  ###### SET INITIAL STATE AND COVARIANCE VALUES #######
  if(is.null(initial.state)){
    stop("Please set provide an initial value for the state and covariance as a named list with entries 'x0' and 'p0'")
  }
  if(!is.list(initial.state)){
    stop("The initial state argument must be a named list with entries 'x0' and 'p0'")
  }
  if(length(initial.state) != 2){
    stop("The initial state should only contain an entry x0 for the initial state, and p0 for the initial covariance")
  }
  if(length(initial.state$x0) != private$number.of.states){
    stop("The length of the initial state for x0 should be equal to the number of states")
  }
  if(any(dim(initial.state$p0) != rep(private$number.of.states,2))){
    stop("The length of the initial state for x0 should be equal to the number of states")
  }
  
  private$set_pred_initial_state(initial.state[["x0"]], initial.state[["p0"]])
  
  return(invisible(self))
}


was_any_data_provided = function(data, self, private)
{
  ###### CHECK DATA #######
  if(is.null(data)){
    if(is.null(private$data)){
      stop("Please supply data.")
    } else {
      message("No new data was provided, reusing the lastest provided data.")
      data = private$data
    }
  }
  
  return(invisible(self))
}

check_and_set_data = function(data, self, private) {
  
  # is data provided, or does private$data hold any data?
  was_any_data_provided(data, self, private)
  
  # convert to data.frame
  data = as.data.frame(data)
  
  # calculate observations in 'call' form e.g. if complex obs lhs e.g:
  # log(y) ~ ... , then we are given y in the data, and must compute log(y)
  data = calculate_complex_observation_lefthandsides(data, self, private)
  
  # Check that inputs, and observations are there
  basic_data_check(data, self, private)
  
  # set timestep
  set_ode_timestep(data, self, private)
  set_simulation_timestep(data, self, private)
  
  # various calculations for tmb's laplace method
  set_data_for_laplace_method(data, self, private)
  
  # only store the obs.names, not the parsed data 
  # example: if we have obs eq log(y) ~ y with name log_y, then we store log_y, but not y itself.
  private$data = data[c(private$obs.names, private$input.names)]
  
  # Return
  return(invisible(self))
}
