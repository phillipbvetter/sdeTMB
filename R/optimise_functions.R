# Functions for construction the AD function from TMB, and calling the
# stats::nlminb estimator for minimizing the negative loglikelihood.

#######################################################
# MAIN CONSTRUCT MAKEADFUN FUNCTION THAT CALL OTHERS
#######################################################

construct_makeADFun = function(self, private){
  
  if(private$method=="laplace"){
    construct_rtmb_laplace_makeADFun(self, private)
  } else {
    construct_kalman_makeADFun(self, private)
  }
  
  return(invisible(self))
}

#######################################################
# CONSTRUCT KALMAN MAKEADFUN
#######################################################
construct_kalman_makeADFun = function(self, private){
  
  ################################################
  # Data
  ################################################
  
  estimation_method = switch(private$method,
                             ekf = 1,
                             ukf = 2
  )
  
  # add mandatory entries to data
  tmb.data = list(
    
    # methods and purpose
    estimation_method = estimation_method,
    ode_solver = private$ode.solver,
    
    # initial
    stateVec = private$initial.state$x0,
    covMat = private$initial.state$p0,
    
    # time-steps
    ode_timestep_size = private$ode.timestep.size,
    ode_timesteps = private$ode.timesteps,
    ode_cumsum_timesteps = private$ode.timesteps.cumsum,
    
    # loss function
    loss_function = private$loss$loss,
    loss_threshold_value = private$loss$c,
    tukey_loss_parameters = private$tukey.pars,
    
    # system size
    number_of_state_eqs = private$number.of.states,
    number_of_obs_eqs = private$number.of.observations,
    number_of_diffusions = private$number.of.diffusions,
    
    # inputs
    inputMat = as.matrix(private$data[private$input.names]),
    
    # observations
    obsMat = as.matrix(private$data[private$obs.names])
  )
  
  # unscented parameters
  ukf_hyperpars_list = list()
  if(private$method=="ukf")
  {
    ukf_hyperpars_list = list(
      ukf_alpha = private$ukf_alpha,
      ukf_beta = private$ukf_beta,
      ukf_kappa = private$ukf_kappa
    )
  }
  
  # MAP Estimation?
  tmb.map.data = list(
    MAP_bool = 0L
  )
  if (!is.null(private$map)) {
    tmb.map.data = list(
      MAP_bool = 1L,
      map_mean__ = private$map$mean[!bool],
      map_cov__ = private$map$cov[!bool,!bool],
      map_ints__ = as.numeric(!bool),
      sum_map_ints__ = sum(as.numeric(!bool))
    )
  }
  
  # construct final data list
  data = c(tmb.data, private$iobs, tmb.map.data, ukf_hyperpars_list)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = lapply(private$parameters, function(x) x[["initial"]]) # Initial parameter values
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = TMB::MakeADFun(data = data,
                       parameters = parameters,
                       map = private$fixed.pars,
                       DLL = private$modelname,
                       silent = TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
}

#######################################################
# CONSTRUCT LAPLACE MAKEADFUN
#######################################################

construct_rtmb_laplace_makeADFun = function(self, private)
{
  
  ################################################
  # Data
  ################################################
  
  # time-steps
  ode_timestep_size = private$ode.timestep.size
  ode_timesteps = private$ode.timesteps
  ode_cumsum_timesteps = private$ode.timesteps.cumsum
  
  # system size
  number_of_state_eqs = private$number.of.states
  number_of_obs_eqs = private$number.of.observations
  number_of_diffusions = private$number.of.diffusions
  
  # inputs
  inputMat = as.matrix(private$data[private$input.names])
  
  # observations
  obsMat = as.matrix(private$data[private$obs.names])
  
  # indices in state parameter vectors corresponding to indices in observations / inputs
  # add 1 because too lazy to change private$iobs from 0-index to 1-indexed.
  iobs = lapply(private$iobs,function(x) x+1)
  
  ################################################
  # Parameters
  ################################################
  
  parameters = c(
    lapply(private$parameters, function(x) x[["initial"]]),
    private$tmb.initial.state.for.parameters
  )
  
  ################################################
  # Define function
  ################################################
  
  elements = get_rtmb_function_elements(self, private)
  
  sde.functions.txt = paste('###### DEFINE FUNCTIONS #######
# drift function
f__ = function(stateVec, parVec, inputVec){
ans = c(F_ELEMENTS)
return(ans)
}

# diffusion function
g__ = function(stateVec, parVec, inputVec){
ans = RTMB::matrix(c(G_ELEMENTS), nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_DIFFUSIONS, byrow=T)
return(ans)
}

h__ = function(stateVec, parVec, inputVec){
ans = c(H_ELEMENTS)
return(ans)
}

hvar__ = function(stateVec, parVec, inputVec){
ans = c(HVAR_ELEMENTS)
return(ans)
}')
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_STATES", 
                                               replacement=deparse(as.numeric(private$number.of.states)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_DIFFUSIONS", 
                                               replacement=deparse(as.numeric(private$number.of.diffusions)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="F_ELEMENTS", 
                                               replacement=paste(elements$f,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="G_ELEMENTS", 
                                               replacement=paste(elements$g,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="H_ELEMENTS", 
                                               replacement=paste(elements$h,collapse=","))
  sde.functions.txt = stringr::str_replace_all(sde.functions.txt, 
                                               pattern="HVAR_ELEMENTS", 
                                               replacement=paste(elements$hvar,collapse=","))
  eval(parse(text=sde.functions.txt))
  
  main.function.txt = paste('laplace.nll = function(p){
# set negative log-likelihood
nll = 0

# small identity matrix
small_identity = diag(1e-8, nrow=NUMBER_OF_STATES, ncol=NUMBER_OF_STATES)

# extract state random effects
stateMat = cbind(STATE_NAMES)

# fixed effects parameter vector
parVec = cbind(FIXED_PARAMETERS)

###### TIME LOOP START #######
for(i in 1:(nrow(obsMat)-1)) {

# Define inputs and use first order input interpolation
inputVec = inputMat[i,]
dinputVec = (inputMat[i+1,] - inputMat[i,])/ode_timesteps[i]

###### BETWEEN TIME POINTS LOOP START #######
for(j in 1:ode_timesteps[i]){

# grab current and next state
x_now = stateMat[ode_cumsum_timesteps[i]+j,]
x_next = stateMat[ode_cumsum_timesteps[i]+j+1,]

# compute drift (vector) and diffusion (matrix)
f = f__(x_now, parVec, inputVec)
g = g__(x_now, parVec, inputVec)
inputVec = inputVec + dinputVec

# assume multivariate gauss distribution according to euler-step
# and calculate the likelihood
z = x_next - (x_now + f * ode_timestep_size[i])
v = (g %*% t(g) + small_identity) * ode_timestep_size[i]
nll = nll - RTMB::dmvnorm(z, Sigma=v, log=TRUE)
}
###### BETWEEN TIME POINTS LOOP END #######
}
###### TIME LOOP END #######

obsMat = RTMB::OBS(obsMat)
###### DATA UPDATE START #######
for(i in 1:NUMBER_OF_OBSERVATIONS){
for(j in 1:length(iobs[[i]])){

# Get index where observation is available
k = iobs[[i]][j]

# Get corresponding input, state and observation
inputVec = inputMat[k,]
stateVec = stateMat[ode_cumsum_timesteps[k]+1,]
obsScalar = obsMat[k,i]

# Observation equation and varianace
h_x = h__(stateVec, parVec, inputVec)[i]
hvar_x = hvar__(stateVec, parVec, inputVec)[i]

# likelihood contribution
nll = nll - RTMB::dnorm(obsScalar, mean=h_x, sd=sqrt(hvar_x), log=TRUE)
}}
###### DATA UPDATE END #######

# return
return(invisible(nll))
}')
  
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_STATES", 
                                               replacement=deparse(as.numeric(private$number.of.states)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_DIFFUSIONS", 
                                               replacement=deparse(as.numeric(private$number.of.diffusions)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="NUMBER_OF_OBSERVATIONS", 
                                               replacement=deparse(as.numeric(private$number.of.observations)))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="STATE_NAMES", 
                                               replacement=paste("p$",private$state.names,sep="",collapse=","))
  main.function.txt = stringr::str_replace_all(main.function.txt, 
                                               pattern="FIXED_PARAMETERS", 
                                               replacement=paste("p$",private$parameter.names,sep="",collapse=","))
  eval(parse(text=main.function.txt))
  
  ################################################
  # Construct Neg. Log-Likelihood
  ################################################
  
  nll = RTMB::MakeADFun(func = laplace.nll, 
                        parameters=parameters, 
                        random=private$state.names,
                        map=private$fixed.pars,
                        silent=TRUE)
  
  # save objective function
  private$nll = nll
  
  # return
  return(invisible(self))
  
}

#######################################################
# OPTIMISE AD FUN
#######################################################

optimize_negative_loglikelihood = function(self, private) {
  
  lower.parameter.bound = unlist(lapply(private$free.pars, function(par) par$lower))
  upper.parameter.bound = unlist(lapply(private$free.pars, function(par) par$upper))
  # If unconstrained optimization was requested remove these bounds
  if(private$unconstrained.optim){
    lower.parameter.bound = -Inf
    upper.parameter.bound = Inf
  }
  
  # IF METHOD IS KALMAN FILTER
  if (any(private$method==c("ekf","ukf"))) {
    
    # use function, gradient and hessian
    if (private$use.hessian) {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           hessian = private$nll$he,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
      # or just function and gradient
    } else {
      comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                           objective = private$nll$fn,
                                                                           gradient = private$nll$gr,
                                                                           lower = lower.parameter.bound,
                                                                           upper = upper.parameter.bound,
                                                                           control=private$control.nlminb))
      )
    }
    
  }
  
  # IF METHOD IS TMB
  if (private$method =="laplace") {
    comptime = system.time( opt <- try_withWarningRecovery(stats::nlminb(start = private$nll$par,
                                                                         objective = private$nll$fn,
                                                                         gradient = private$nll$gr,
                                                                         lower = lower.parameter.bound,
                                                                         upper = upper.parameter.bound,
                                                                         control=private$control.nlminb))
    )
  }
  
  # DID THE OPTIMIZATION FAIL?
  if (inherits(opt,"try-error")) {
    
    message("The optimisation failed due to the following error: \n\n\t",opt)
    
    if(stringr::str_detect(opt,"NA/NaN")){
      
      message("You should consider the following to circumvent the error:
              1. Reduce the ODE step-size via argument 'ode.timestep'.
              2. Run the optimization with / without the hessian via argument 'use.hessian'.
              3. Change / explore parameter initial values.
              4. Extract the function handlers with the 'construct_nll' method and investigate outputs, or try other optimizers.
              5. Change the optimization tolerances for 'nlminb' with the 'control' argument.")
      
    }
    
    private$opt = NULL
    
    return(invisible(self))
  }
  
  # store optimization object
  private$opt = opt
  
  # extract maxmimum gradient component, and format computation time to 5 digits
  outer_mgc = max(abs(private$nll$gr(opt$par)))
  comp.time = format(round(as.numeric(comptime["elapsed"])*1e4)/1e4,digits=5,scientific=F)
  
  # print convergence and timing result
  if(outer_mgc > 1){
    message("BEWARE: THE MAXIMUM GRADIENT COMPONENT APPEARS TO BE LARGE ( > 1 ) - THE FOUND OPTIMUM MIGHT BE INVALID.")
  }
  message("\t Optimization finished!:
            Elapsed time: ", comp.time, " seconds.
            The objective value is: ",format(opt$objective,scientific=T),"
            The maximum gradient component is: ",format(outer_mgc,digits=2,scientific=T),"
            The convergence message is: ", opt$message,"
            Iterations: ",opt$iterations,"
            Evaluations: Fun: ",opt$evaluations["function"]," Grad: ",opt$evaluations[["gradient"]],"
            See stats::nlminb for available tolerance/control arguments."
  )
  
  # For TMB method: run sdreport
  if (private$method=="laplace") {
    message("Calculating random effects standard deviation...")
    private$sdr = RTMB::sdreport(private$nll)
    message("Finished!")
  }
  
  # return
  return(invisible(self))
}

#######################################################
# MAKE RETURN DATA NICE AFTER OPTIMIZATION
#######################################################

create_return_fit = function(self, private, calculate.laplace.onestep.residuals) {
  
  if (is.null(private$opt)) {
    return(NULL)
  }
  
  # store the provided data in the fit
  private$fit$data = private$data
  
  # get convergence
  private$fit$convergence = private$opt$convergence
  
  # store the object in fit - these gives access to e.g.
  private$fit$.__object__ = self$clone()
  
  ################################################
  # FOR KALMAN FILTERS
  ################################################
  
  if (any(private$method == c("ekf","ukf"))) {
    
    
    ################################################
    # BASICS
    ################################################
    
    # objective value
    private$fit$nll = private$opt$objective
    
    # gradient
    private$fit$nll.gradient = try_withWarningRecovery(
      {
        nll.grad = as.vector(private$nll$gr(private$opt$par))
        names(nll.grad) = names(private$free.pars)
        nll.grad
      }
    )
    if (inherits(private$fit$nll.gradient,"try-error")) {
      private$fit$nll.gradient = NA
    }
    
    # hessian
    private$fit$nll.hessian = try_withWarningRecovery(
      {
        nll.hess = private$nll$he(private$opt$par)
        rownames(nll.hess) = names(private$free.pars)
        colnames(nll.hess) = names(private$free.pars)
        nll.hess
      }
    )
    if (inherits(private$fit$nll.hessian, "try-error")) {
      private$fit$nll.hessian = NA
    }
    
    # parameter estimates and standard deviation
    private$fit$par.fixed = private$opt$par
    private$fit$sd.fixed = tryCatch(sqrt(diag(solve(private$nll$he(private$opt$par)))),
                                    error=function(e) NA,
                                    warning=function(w) NA
    )
    
    # parameter covariance matrix by solving the inverse hessian
    private$fit$cov.fixed = tryCatch(solve(private$nll$he(private$opt$par)),
                                     error=function(e) NA,
                                     warning=function(w) NA
    )
    
    ################################################
    # STATES, RESIDUALS, OBSERVATIONS ETC.
    ################################################
    
    # Extract reported items from nll
    rep = private$nll$report()
    
    # Prior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, do.call(rbind,rep$xPrior)))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPrior,diag)))))
    
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$prior = as.data.frame(temp.states)
    private$fit$states$sd$prior = as.data.frame(temp.sd)
    private$fit$states$cov$prior = rep$pPrior
    names(private$fit$states$cov$prior) = paste("t = ",private$data$t,sep="")
    
    # Posterior States
    temp.states = try_withWarningRecovery(cbind(private$data$t, do.call(rbind,rep$xPost)))
    temp.sd = try_withWarningRecovery(cbind(private$data$t, sqrt(do.call(rbind,lapply(rep$pPost,diag)))))
    colnames(temp.states) = c("t",private$state.names)
    colnames(temp.sd) = c("t",private$state.names)
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)
    private$fit$states$cov$posterior = rep$pPost
    names(private$fit$states$cov$posterior) = paste("t = ",private$data$t,sep="")
    
    # Residual
    # rowNAs = as.matrix(!is.na(do.call(cbind, private$data[private$obs.names]))[-1,])
    rowNAs = as.matrix(!is.na(private$data[private$obs.names])[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    innovation = rep$Innovation
    innovation.cov = rep$InnovationCovariance
    innovation[[1]] = NULL
    innovation.cov[[1]] = NULL
    
    temp.res = matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    temp.var =  matrix(nrow=length(private$data$t)-1, ncol=private$number.of.observations)
    
    # do.call(rbind, lapply(rep$Innovation, "length<-", private$m))
    for (i in seq_along(private$data$t[-1])) {
      if (sumrowNAs[i] > 0) {
        temp.res[i,rowNAs[i,]] = innovation[[i]]
        temp.var[i,rowNAs[i,]] = diag(innovation.cov[[i]])
      }
    }
    temp.res = cbind(private$data$t[-1], temp.res)
    temp.sd = cbind(private$data$t[-1], sqrt(temp.var))
    
    names(innovation.cov) = paste("t = ",private$data$t[-1],sep="")
    
    # should we remove the empty matrices?
    # innovation.cov = innovation.cov[sumrowNAs!=0]
    
    colnames(temp.res) = c("t",private$obs.names)
    colnames(temp.sd) = c("t",private$obs.names)
    private$fit$residuals$mean = as.data.frame(temp.res)
    private$fit$residuals$sd = as.data.frame(temp.sd)
    private$fit$residuals$normalized = as.data.frame(temp.res)
    private$fit$residuals$normalized[,-1] = private$fit$residuals$normalized[,-1]/temp.sd[,-1]
    private$fit$residuals$cov = innovation.cov
    
    
    # Observations
    # We need all states, inputs and parameter values to evaluate the observation
    # put them in a list
    listofvariables.prior = c(
      as.list(private$fit$states$mean$prior[-1]),
      as.list(private$fit$par.fixed),
      as.list(private$fit$data)
    )
    listofvariables.posterior = c(
      as.list(private$fit$states$mean$posterior[-1]),
      as.list(private$fit$par.fixed),
      as.list(private$fit$data)
    )
    obs.df.prior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.prior)})
    )
    obs.df.posterior = as.data.frame(
      lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = listofvariables.posterior)})
    )
    private$fit$observations$mean$prior = data.frame(t=private$data$t, obs.df.prior)
    private$fit$observations$mean$posterior = data.frame(t=private$data$t, obs.df.posterior)
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  }
  
  ################################################
  # FOR TMB
  ################################################
  
  if (private$method == "laplace") {
    
    n = private$number.of.states
    
    # Objective and Gradient
    private$fit$nll = private$opt$objective
    private$fit$nll.gradient = private$sdr$gradient.fixed
    names(private$fit$nll.gradient) = names(private$free.pars)
    
    # Parameter Estimate
    private$fit$par.fixed = private$opt$par
    private$fit$sd.fixed = diag(private$sdr$cov.fixed)
    
    # Parameter Covariance
    private$fit$cov.fixed = private$sdr$cov.fixed
    
    # Posterior States (Smoothed)
    temp = cbind(matrix(private$sdr$par.random, ncol=n),
                 matrix(sqrt(private$sdr$diag.cov.random), ncol=n))[private$ode.timesteps.cumsum+1, ]
    temp.states = cbind(private$data$t, matrix(temp[,1:n],nrow=length(private$data$t)))
    temp.sd = cbind(private$data$t, matrix(temp[,(n+1):(2*n)],nrow=length(private$data$t)))
    #
    private$fit$states$mean$posterior = as.data.frame(temp.states)
    private$fit$states$sd$posterior = as.data.frame(temp.sd)
    colnames(private$fit$states$sd$posterior) = c("t",private$state.names)
    colnames(private$fit$states$mean$posterior) = c("t",private$state.names)
    
    # Residuals
    rowNAs = as.matrix(!is.na(do.call(cbind,private$data[private$obs.names]))[-1,])
    sumrowNAs = rowSums(rowNAs)
    
    # compute one-step residuals
    if(calculate.laplace.onestep.residuals){
      message("Calculating one-step ahead residuls...")
      osa = RTMB::oneStepPredict(private$nll,
                                 observation.name="obsMat",
                                 method="oneStepGaussian",
                                 trace=TRUE)
      private$fit$one_step_residuals = osa
    }
    
    # t-values and Pr( t > t_test )
    private$fit$tvalue = private$fit$par.fixed / private$fit$sd.fixed
    private$fit$Pr.tvalue = 2*pt(q=abs(private$fit$tvalue),df=sum(sumrowNAs),lower.tail=FALSE)
    
  }
  
  # Set S3 class
  class(private$fit) = "sdem.fit"
  
  return(invisible(self))
}

#######################################################
# CREATE DATA.FRAME FOR PREDICT METHOD
#######################################################

# construct_predict_dataframe = function(pars, rep, data, return.covariance, return.k.ahead, self, private){
#   
#   # Simlify variable names
#   n = private$number.of.states
#   n.ahead = private$n.ahead
#   state.names = private$state.names
#   last.pred.index = private$last.pred.index
#   
#   df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
#   disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
#   disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
#   var_bool = !stringr::str_detect(disp_names,"cor")
#   if(return.covariance){
#     disp_names = sprintf(rep("cov.%s.%s",n^2),rep(state.names,each=n),rep(state.names,n))
#     disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
#   }
#   names(df.out) = c("i.","j.","t.i","t.j","k.ahead",state.names,disp_names)
#   ran = 0:(last.pred.index-1)
#   df.out["i."] = rep(ran,each=n.ahead+1)
#   df.out["j."] = df.out["i."] + rep(0:n.ahead,last.pred.index)
#   df.out["t.i"] = rep(data$t[ran+1],each=n.ahead+1)
#   df.out["t.j"] = data$t[df.out[,"i."]+1+rep(0:n.ahead,last.pred.index)]
#   df.out["k.ahead"] = rep(0:n.ahead,last.pred.index)
#   df.out[,state.names] = do.call(rbind,rep$xk__)
#   if(return.covariance){
#     df.out[,disp_names] = do.call(rbind,rep$pk__)
#   } else {
#     df.out[,disp_names] = do.call(rbind, lapply(rep$pk__,function(x) do.call(rbind, apply(x,1, function(y) as.vector(cov2cor(matrix(y,ncol=2,byrow=T))),simplify=FALSE))))
#     diag.ids = seq(from=1,to=n^2,by=n+1)
#     df.out[,disp_names[diag.ids]] = do.call(rbind,rep$pk__)[,diag.ids]
#   }
#   
#   # calculate observations at every time-step in predict
#   inputs.df = private$data[df.out[,"j."]+1,private$input.names]
#   
#   env.list = c(
#     as.list(df.out[state.names]),
#     as.list(inputs.df),
#     as.list(pars)
#   )
#   
#   obs.df.predict = as.data.frame(
#     lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = env.list)})
#   )
#   names(obs.df.predict) = paste(private$obs.names,".predict",sep="")
#   df.out = cbind(df.out, obs.df.predict)
#   
#   # add data observation to output data.frame 
#   obs.df = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
#   names(obs.df) = paste(private$obs.names,".data",sep="")
#   df.out = cbind(df.out, obs.df)
#   
#   # return only specific n.ahead
#   df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
#   class(df.out) = c("sdem.pred", "data.frame")
#   
#   return(df.out)
# }

construct_predict_rcpp_dataframe = function(pars, predict.list, data, return.covariance, return.k.ahead, self, private){
  
  # Simlify variable names
  n = private$number.of.states
  n.ahead = private$n.ahead
  state.names = private$state.names
  last.pred.index = private$last.pred.index
  
  # Create return data.frame
  df.out = data.frame(matrix(nrow=last.pred.index*(n.ahead+1), ncol=5+n+n^2))
  disp_names = sprintf(rep("cor.%s.%s",n^2),rep(state.names, each=n),rep(state.names,n))
  disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  var_bool = !stringr::str_detect(disp_names,"cor")
  if(return.covariance){
    disp_names = sprintf(rep("cov.%s.%s",n^2),rep(state.names,each=n),rep(state.names,n))
    disp_names[seq.int(1,n^2,by=n+1)] = sprintf(rep("var.%s",n),state.names)
  }
  names(df.out) = c("i.","j.","t.i","t.j","k.ahead",state.names,disp_names)
  
  # Fill out data.frame
  ran = 0:(last.pred.index-1)
  df.out["i."] = rep(ran,each=n.ahead+1)
  df.out["j."] = df.out["i."] + rep(0:n.ahead,last.pred.index)
  df.out["t.i"] = rep(data$t[ran+1],each=n.ahead+1)
  df.out["t.j"] = data$t[df.out[,"i."]+1+rep(0:n.ahead,last.pred.index)]
  df.out["k.ahead"] = rep(0:n.ahead,last.pred.index)
  df.obs = df.out[c("i.","j.","t.i","t.j","k.ahead")]
  
  ##### STATES PREDICTIONS ######
  df.out[,state.names] = do.call(rbind,lapply(predict.list$Xpred, function(cur.list) do.call(rbind, cur.list)))
  if(return.covariance){ #covariance
    df.out[,disp_names] = do.call(rbind, lapply(predict.list$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(x))) } ))
  } else { #correlation
    df.out[,disp_names] = do.call(rbind, lapply(predict.list$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(cov2cor(x)))) } )) 
    diag.ids = seq(from=1,to=n^2,by=n+1)
    df.out[,disp_names[diag.ids]] = do.call(rbind, lapply(predict.list$Ppred, function(cur.list){ do.call(rbind, lapply(cur.list, function(x) as.vector(diag(x)))) } ))
  }
  
  ##### OBSERVATION PREDICTIONS #####
  # calculate observations at every time-step in predict
  inputs.df = private$data[df.out[,"j."]+1,private$input.names]
  
  # create environment
  env.list = c(
    as.list(df.out[state.names]),
    as.list(inputs.df),
    as.list(pars)
  )
  
  # calculate observations
  obs.df.predict = as.data.frame(
    lapply(private$obs.eqs.trans, function(ls){eval(ls$rhs, envir = env.list)})
  )
  names(obs.df.predict) = paste(private$obs.names,".predict",sep="")
  # df.out = cbind(df.out, obs.df.predict)
  
  # add data observation to output data.frame 
  obs.df.data = private$data[df.out[,"j."]+1, private$obs.names, drop=F]
  names(obs.df.data) = paste(private$obs.names,".data",sep="")
  # df.out = cbind(df.out, obs.df)
  
  df.obs = cbind(df.obs, obs.df.predict, obs.df.data)
  
  # return only specific n.ahead
  df.out = df.out[df.out[,"k.ahead"] %in% return.k.ahead,]
  df.obs = df.obs[df.obs[,"k.ahead"] %in% return.k.ahead,]
  
  list.out = list(states = df.out, observations = df.obs)
  class(list.out) = c(class(list.out), "sdem.pred")
  
  return(list.out)
}

construct_simulate_rcpp_dataframe = function(pars, predict.list, data, return.covariance, return.k.ahead, self, private){
  
  
  list.out = vector("list",length=private$number.of.states)
  names(list.out) = private$state.names
  
  setRownames = function(obj, nm){rownames(obj) = nm; return(obj)}
  setColnames = function(obj, nm){colnames(obj) = nm; return(obj)}
  
  
  for(i in seq_along(list.out)){
    list.out[[i]] = stats::setNames(
      lapply(predict.list, function(ls.outer){
        # setRownames(
          t(do.call(cbind, lapply(ls.outer, function(ls.inner) ls.inner[,i])))
          # paste0("k.ahead", 0:private$n.ahead)
        # )
      }),
      paste0("t", head(data$t, private$last.pred.index))
    )
  }
  
  # Compute the prediction times for each horizon
  ran = 0:(private$last.pred.index-1)
  t.j = data$t[rep(ran,each=private$n.ahead+1)+1+rep(0:private$n.ahead,private$last.pred.index)]
  t.j.splitlist = split(t.j, ceiling(seq_along(t.j)/(private$n.ahead+1)))
  list.of.time.vectors = lapply(t.j.splitlist, function(x) data.frame(t.j=x))
  names(list.of.time.vectors) = names(list.out[[1]])
  
  list.out2 = c(list.out, list(prediction_times = list.of.time.vectors))
  
  return(list.out2)
}


