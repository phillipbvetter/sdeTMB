#' @title R6 Class constructor for sdemTMB
#' @description Call sdemTMB$new() to construct a new instance of the class object
#' for specifying a new stochastic state space model.
#' @returns A model object of class 'sdemTMB'
#' @export
sdemTMB = R6::R6Class(
  classname = "sdemTMB",
  # public fields
  public = list(

    #' @description Construct new sdemTMB model object
    initialize = function() {
      message("Created New SDEM Object")
      # public fields
      #
      # modelname
      private$modelname = "sde_model"
      private$cppfile.directory = getwd()
      private$cppfile.path = paste(getwd(),"/",private$modelname,sep="")
      # model equations
      private$sys.eqs = NULL
      private$obs.eqs = NULL
      private$obs.var = NULL
      private$alg.eqs = NULL
      private$inputs = list(list(input=quote(t)))
      names(private$inputs) = "t"
      private$constants = NULL
      private$parameters = NULL
      private$initial.state = NULL
      private$tmb.initial.state.for.parameters = NULL
      private$iobs = NULL
      # after algebraics
      private$sys.eqs.trans = NULL
      private$obs.eqs.trans = NULL
      private$obs.var.trans = NULL
      # options
      private$method = "ekf"
      private$use.hessian = FALSE
      private$state.dep.diff = FALSE
      private$lamperti = list(transform="identity",states=NULL)
      private$compile = FALSE
      private$loss = list(loss=0L,c=3)
      private$tukey.pars = rep(0,4)
      private$silent = FALSE
      private$map = NULL
      # hidden
      private$linear = NULL
      private$fixed.pars = NULL
      private$free.pars = NULL
      private$trigger = NULL
      # names
      private$state.names = NULL
      private$obs.names = NULL
      private$obsvar.names = NULL
      private$input.names = "t"
      private$parameter.names = NULL
      private$constant.names = NULL
      # lengths
      private$n = 0
      private$m = 0
      private$ng = 0
      # differentials
      private$diff.processes = NULL
      private$diff.terms = NULL
      # build flag
      private$build = FALSE
      # data, nll, opt
      private$data = NULL
      private$nll = NULL
      private$opt = NULL
      private$fit = NULL
    },

    #' @description Add system equations to model object.
    #'
    #' Example:
    #'
    #' \code{form = dx ~ theta * (mu - x) * dt + sigma * dw}.
    #' @param form formula class specifying the stochastic differential equation(s) to be added to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple SDEs at the same time.
    add_systems = function(form,...) {
      # lapply for multi-arguments
      lapply(c(form,...), function(form) {
        # Function Body
        res = check_system_eqs(form, self, private)
        check_name(names(res), "state", self, private)
        private$sys.eqs[[names(res)]] = res[[1]]
        private$state.names = names(private$sys.eqs)
        # Function Body
      })
      # private$n = length(private$sys.eqs)
      # if model was already built, set compile flag to true
      was_model_already_built(self, private)
      #
      return(invisible(self))
    },

    #' @description Add observation equations to model object.
    #'
    #' Example:
    #'
    #' Assuming that \code{x} is a state defined through \code{add_systems} with \code{form = dx ~ ...} then define an observation related to
    #' \code{x} as \code{y ~ x}.
    #' @param form formula class specifying the obsevation equation to be added to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple observation equations at a time.
    add_observations = function(form,...) {
      # lapply for multi-arguments
      lapply(c(form,...), function(form) {
        # Function Body
        res = check_obsservation_eqs(form, self, private)
        check_name(names(res), "obs", self, private)
        private$obs.eqs[[names(res)]] = res[[1]]
        private$obs.names = names(private$obs.eqs)
        # Function Body
      })
      # if model was already built, set compile flag to true
      was_model_already_built(self, private)
      #
      return(invisible(self))
    }

  ),
  # private fields
  private = list(
    # modelname
    modelname = character(0),
    cppfile.directory = NULL,
    cppfile.path = character(0),
    # model equations
    sys.eqs = NULL,
    obs.eqs = NULL,
    obs.var = NULL,
    alg.eqs = NULL,
    inputs = NULL,
    constants = NULL,
    parameters = NULL,
    initial.state = NULL,
    tmb.initial.state.for.parameters = NULL,
    iobs = NULL,
    # after algebraics
    sys.eqs.trans = NULL,
    obs.eqs.trans = NULL,
    obs.var.trans = NULL,
    # names
    state.names = NULL,
    obs.names = NULL,
    obsvar.names = NULL,
    input.names = NULL,
    parameter.names = NULL,
    constant.names = NULL,
    # options
    method = NULL,
    use.hessian = NULL,
    state.dep.diff = NULL,
    lamperti = NULL,
    compile = NULL,
    loss = NULL,
    tukey.pars = NULL,
    silent = NULL,
    map = NULL,
    ode.timestep = NULL,
    ode.dt = NULL,
    ode.N = NULL,
    ode.Ncumsum = NULL,
    # hidden
    linear = NULL,
    fixed.pars = NULL,
    free.pars = NULL,
    trigger = NULL,
    # lengths
    n = NULL,
    m = NULL,
    ng = NULL,
    # differentials
    diff.processes = NULL,
    diff.terms = NULL,
    # build flag
    build = NULL,
    # data, nll, opt
    data = NULL,
    nll = NULL,
    opt = NULL,
    sdr = NULL,
    fit = NULL
  )
)

sdemTMB$set("public","add_trans_systems",
            function(form) {
              res = check_system_eqs(form, self, private, silent=T)
              private$sys.eqs.trans[[names(res)]] = res[[1]]
              return(invisible(self))
            }
)

sdemTMB$set("public","add_trans_observations",
            function(form) {
              res = check_obsservation_eqs(form, self, private, silent=T)
              private$obs.eqs.trans[[names(res)]] = res[[1]]
              return(invisible(self))
            }
)

sdemTMB$set("public","add_observation_variances",
            function(form,...) {
              # lapply for multi-arguments
              lapply(c(form,...), function(form) {
                # Function Body
                res = check_observation_variance_eqs(form, self, private)
                check_name(names(res), "obsvar", self, private)
                private$obs.var[[names(res)]] = res[[1]]
                private$obsvar.names = names(private$obs.var)
                # Function Body
              })
              # if model was already built, set compile flag to true
              was_model_already_built(self, private)
              #
              return(invisible(self))
            }
)

sdemTMB$set("public","add_trans_observation_variances",
            function(form) {
              res = check_observation_variance_eqs(form, self, private, silent=T)
              private$obs.var.trans[[names(res)]] = res[[1]]
              return(invisible(self))
            }
)

sdemTMB$set("public","add_inputs",
            function(...) {
              args = as.list(match.call()[-1])
              # lapply for multi-arguments
              lapply(args, function(args) {
                # Function Body
                res = check_inputs(args, self, private)
                check_name(names(res), "input", self, private)
                private$trigger = is_this_a_new_name(names(res),private$input.names)
                private$inputs[[names(res)]] = res[[1]]
                private$input.names = names(private$inputs)
                # Function Body
              })
              # if model was already built, set compile flag to true
              if (private$trigger) {
                was_model_already_built(self, private)
              }
              #
              return(invisible(self))
            }
)

sdemTMB$set("public","add_parameters",
            function(...,parameter.matrix=NULL) {

              # single parameters
              parameter.forms = list(...)
              if(length(parameter.forms)>0){
                lapply(parameter.forms, function(parameter.form) {
                  check_parameter(parameter.form)
                  parname = deparse1(parameter.form[[2]])
                  par.val.lb.ub = eval(parameter.form[[3]])
                  check_name(parname, "pars", self, private)
                  private$trigger = is_this_a_new_name(parname, private$parameter.names)
                  private$parameters[[parname]] = list(init=par.val.lb.ub[1],
                                                       lb=par.val.lb.ub[2],
                                                       ub=par.val.lb.ub[3])
                  # fixed parameters
                  if(is.na(par.val.lb.ub[2]) & is.na(par.val.lb.ub[3])){
                    private$fixed.pars[[parname]] = factor(NA)
                  }
                })
                private$parameter.names = names(private$parameters)
              }

              # parameter matrix
              if(!is.null(parameter.matrix)){
                check_parameter_matrix(parameter.matrix, self, private)
                parnames = rownames(parameter.matrix)
                lapply( 1:nrow(parameter.matrix), function(i) {
                  parname = parnames[i]
                  par.val.lb.ub = parameter.matrix[i,]
                  check_name(parname, "pars", self, private)
                  private$trigger = is_this_a_new_name(parname, private$parameter.names)
                  private$parameters[[parname]] = list(init=par.val.lb.ub[1],
                                                       lb=par.val.lb.ub[2],
                                                       ub=par.val.lb.ub[3])
                  # fixed parameters
                  if(is.na(par.val.lb.ub[2]) & is.na(par.val.lb.ub[3])){
                    private$fixed.pars[[parname]] = factor(NA)
                  }
                })
                private$parameter.names = names(private$parameters)
              }

              # if model was already built set compile flag
              if (private$trigger) {
                was_model_already_built(self, private)
              }

              # return
              return(invisible(self))
            }
)

sdemTMB$set("public","add_parameters2",
            function(parameters) {

              parnames = check_parameter_matrix(parameters, self, private)[[1]]
              lapply(parnames, function(parname) check_name(parname, "pars", self, private))
              private$trigger = any(sapply(parnames, function(parname) is_this_a_new_name(parname,private$parameter.names)))
              for (i in 1:nrow(parameters)) {
                private$parameters[[parnames[i]]] = list(init=parameters[i,1],
                                                         lb=parameters[i,2],
                                                         ub=parameters[i,3])
              }
              private$parameter.names = names(private$parameters)
              # store fixed parameters (those with NA bounds)
              bool = unlist(lapply(private$parameters, function(x) all(is.na(c(x[["lb"]],x[["ub"]])))))
              fixed.pars.names = private$parameter.names[bool]
              for(name in fixed.pars.names) {
                private$fixed.pars[[name]] = factor(NA)
              }

              # if model was already built set compile flag to true
              if (private$trigger) {
                was_model_already_built(self, private)
              }

              # return
              return(invisible(self))
            }
)

sdemTMB$set("public","add_algebraics",
            function(form,...) {
              lapply(c(form,...), function(form) {
                # Function Body
                res = check_algebraics(form, self, private)
                private$alg.eqs[[names(res)]] = res[[1]]
                # Function Body
              })
              # if model was already built, set compile flag to true
              was_model_already_built(self, private)
              #
              return(invisible(self))
            }
)

sdemTMB$set("public","add_constants",
            function(form,...) {
              lapply(c(form,...), function(form) {
                # Function Body
                res = check_constants(form, self, private)
                check_name(names(res), "constants", self, private)
                private$trigger = is_this_a_new_name(names(res),private$constant.names)
                private$constants[[names(res)]] = res[[1]]
                private$constant.names = names(private$constants)
                # Function Body
              })
              # if model was already built, set compile flag to true
              if (private$trigger) {
                was_model_already_built(self, private)
              }
              #
              return(invisible(self))
            }
)

#### ################ ################ ############
#### ################ ################ ############
#### ################ ################ ############

sdemTMB$set("public","set_initial_state",
            function(mean,cov) {
              if (is.null(private$sys.eqs)) {
                stop("Please specify system equations first")
              }
              if (!is.numeric(mean)) {
                stop("The mean vector is not a numeric")
              }
              if (any(is.na(mean))) {
                stop("The mean vector contains NAs.")
              }
              if (length(mean)!=length(private$sys.eqs)) {
                stop("The initial state vector should have length ",length(private$sys.eqs))
              }
              if (!all(dim(cov)==c(length(private$sys.eqs),length(private$sys.eqs)))) {
                stop("The covariance matrix should be square with dimension ", length(private$sys.eqs))
              }
              if (!is.numeric(cov)) {
                stop("The covariance matrix is not a numeric")
              }
              if (any(is.na(cov))) {
                stop("The covariance matrix contains NAs")
              }

              private$initial.state = list(mean=mean,cov=as.matrix(cov))
              return(invisible(self))
            }
)

sdemTMB$set("public","set_method",
            function(method) {
              # is the method a string?
              if (!(is.character(method))) {
                stop("You must pass a string")
              }
              # choose one of the available methods
              available_methods = c("ekf","ekf2","ukf","tmb","tmb_exact")
              if (!(method %in% available_methods)) {
                stop("That method is not available. Please choose one of the following instead: \n
                  1. Extended Kalman Filter - method = 'ekf'
                  2. Unscented Kalman Filter - method = 'ukf'
                  3. Laplace Approx using Random Effects - method = 'tmb',
                  4. Exact using Random Effects (Linear SDEs) - method = 'tmb_exact'")
              }
              private$method = method
              return(invisible(self))
            })

sdemTMB$set("public","set_lamperti",
            function(transform,states=NULL) {
              # must be a string
              if (!(is.character(transform))) {
                stop("You must pass a string")
              }
              if (!is.null(states)){
                if (!(is.character(states))) {
                  stop("You must pass a vector of state names")
                }
                bool = states %in% names(private$sys.eqs)
                if (!all(bool)) {
                  stop("The following state names don't exist: \n\t ",states[!bool])
                }
              }
              # available transforms
              available_transforms = c("identity","log","logit","sqrt-logit","power-logit")
              if (!(transform %in% available_transforms)) {
                stop("That method is not available. Please choose one of the following instead: \n
                  1. For 'dw' use no transform = 'identity' (default)
                  2. For 'x * dw' use transform = 'log'
                  3. For 'x * (1 - x) * dw' use transform = 'logit'
                  4. For 'sqrt( x * (1 - x) ) * dw' use transform = 'sqrt-logit'
                  5. For 'x * (1 - x^a) * dw' for a>0 use transform 'power-logit'")
              }
              private$lamperti = list(transform=transform,states=states)
              return(invisible(self))
            })

sdemTMB$set("public","set_loss",
            function(loss,c=3) {
              # is the method a string?
              if (!(is.character(loss))) {
                stop("You must pass a string")
              }
              # choose one of the available methods
              available_losses = c("standard","huber","tukey")
              if (!(loss %in% available_losses)) {
                stop("That method is not available. Please choose one of the following instead: \n
                  1. Quadratic loss = 'standard'
                  2. Quadratic-Linear loss = 'huber'
                  3. Quadratic-Constant loss = 'tukey'")
              }
              if (loss=="tukey") {
                l = 1L
                # compute tukey approx coefficients
                rtukey = seq(0,100,by=1e-2)
                ctukey = c
                funtukey = function(r){
                  ifelse(r^2 <= ctukey^2,
                         ctukey^2/6 * (1-(1-(r/ctukey)^2)^3),
                         ctukey^2/6
                  )
                }
                tukeyloss = function(.pars){
                  res = sum((funtukey(rtukey) - .pars[4]*(sigmoid(rtukey,a=.pars[1],b=.pars[2])+.pars[3]))^2)
                }
                tukeyopt = nlminb(start=rep(1,4),objective=tukeyloss)
                private$tukey.pars = tukeyopt$par
              } else if (loss=="huber") {
                l = 2L
              } else {
                l = 0L
              }
              private$loss = list(loss=l,c=c)
              return(invisible(self))
            })

sdemTMB$set("public","set_modelname",
            function(name) {
              # was a string passed?
              if (!is.character(name)) {
                stop("The modelname must be a string")
              }
              private$modelname = name
              private$cppfile.path = paste(private$cppfile.directory,"/",name,sep="")
            }
)

sdemTMB$set("public","set_cppfile_directory",
            function(directory) {
              # was a string passed?
              if (!is.character(directory)) {
                stop("You must pass a string")
              }
              # does the exist?
              if (!dir.exists(directory)) {
                stop("The specified directory does not exist")
              }
              private$cppfile.directory = directory

              # update private$cppfile.path by calling set_modelname
              self$set_modelname(private$modelname)
              return(invisible(self))
            })

sdemTMB$set("public","set_compile",
            function(bool) {
              # is bool logical
              if (!is.logical(bool)) {
                stop("You must pass TRUE or FALSE")
              }
              private$compile = bool
              return(invisible(self))
            }
)

sdemTMB$set("public","set_silence",
            function(bool) {
              # is bool logical
              if (!is.logical(bool)) {
                stop("You must pass TRUE or FALSE")
              }
              private$silent = bool
              return(invisible(self))
            }
)

sdemTMB$set("public","set_map",
            function(mean,cov) {
              if (!is.numeric(mean)) {
                stop("The MAP mean vector is not numeric")
              }
              if (length(mean)!=length(private$parameters)) {
                stop("The MAP parameter vector should have length ",length(private$parameters))
              }
              if (!is.matrix(cov)) {
                stop("The MAP covariance matrix is not of class matrix")
              }
              if (!all(dim(cov)==rep(length(private$parameters),2))) {
                stop("The MAP covariance matrix should be square with dimension ", length(private$parameters))
              }
              private$map = list(mean=mean,cov=cov)
            }
)

sdemTMB$set("public","set_timestep",
            function(dt) {
              # OK if NULL
              if(is.null(dt)){
                private$ode.timestep = dt
                return(invisible(NULL))
              }
              # must be numeric
              if (!is.numeric(dt)) {
                stop("The timestep should be a numeric value.")
              }
              # must have length 1
              if (length(dt)!=1) {
                stop("Timestep must have length 1.")
              }
              private$ode.timestep = dt
            }
)


sdemTMB$set("public","use_hessian",
            function(bool) {
              if (!is.logical(bool)) {
                stop("This must be a logical")
              }
              private$use.hessian = bool
            }
)

# sdemTMB$set("public","set_tmb_init_state",
#          function(state_vec) {
#
#            if (any(private$method != c("tmb","tmb_exact"))) {
#              stop("This option is only relevant if you use a tmb random effects style estimation. Use set_method first.")
#            }
#            if (!is.list(state_vec)) {
#              stop("You must pass a named-list of numeric values")
#            }
#            if (any(sapply(state_vec, function(x) !is.numeric(x)))){
#              stop("The list entries must contain numeric values")
#            }
#            if (!all(private$state.names==names(state_vec))) {
#              stop("The names of the passed list do not match the state names")
#            }
#            private$tmb.initial.state = state_vec
#          }
# )

#### ################ ################ ############
#### ################ ################ ############
#### ################ ################ ############

sdemTMB$set("public","build_model",
            function() {

              # basic checks for model, add class n, ng, m, diff procs
              init_build(self, private)

              # apply algebraics
              check_algebraics_before_applying(self, private)
              apply_algebraics(self, private)

              # update diff.terms and apply lamperti
              update_diffterms(self, private)
              apply_lamperti(self, private)
              update_diffterms(self, private)

              # check if model is ok
              lastcheck_before_compile(self, private)

              # compile cpp file
              compile_cppfile(self, private)

              # set build
              private$build = TRUE

              return(invisible(self))
            }
)

sdemTMB$set("public","estimate",
            function(data, return.fit=TRUE, return.nll=FALSE, use.hessian=FALSE,
                     ode.timestep=NULL, silence=FALSE, compile=FALSE) {

              # set flags
              self$use_hessian(use.hessian)
              self$set_timestep(ode.timestep)
              self$set_silence(silence)
              self$set_compile(compile)

              # if the model isnt built we must build
              if (!private$build | private$compile) {
                message("Building model...")
                self$build_model()
              }

              # check and set data
              check_and_set_data(data, self, private)

              # construct neg. log-likelihood function
              optlist = construct_and_optimise(self, private, return.fit, return.nll)

              # if return.nll just return the function objective and exit
              if(return.nll){
                message("Returning AD objective function, and exiting...")
                return(private$nll)
              }

              # create and return fit object
              if(return.fit){
                create_return_fit(self, private)
                return(private$fit)
              }

              # return optimization and cpu-time instead of fit
              if(!return.fit){
                return(optlist)
              }
            }
)


#### ################ ################ ############
#### ################ ################ ############
#### ################ ################ ############

sdemTMB$set("public","print",
            function() {
              n = length(private$sys.eqs)
              m = length(private$obs.eqs)
              p = length(private$inputs)-1
              ng = max(length(unique(unlist(lapply(private$sys.eqs, function(x) x$diff))))-1,0)
              q = length(private$alg.eqs)
              par = length(private$parameters)
              fixedpars = length(private$fixed.pars)
              # If the model is empty
              cat("A stochastic state space model called",private$modelname,"\n")
              basic.data = c(n,ng,m,p,par)
              row.names = c("Number of States","Number of Diffusions",
                            "Number of Observations","Number of Inputs",
                            "Number of Parameters")
              mat=data.frame(basic.data,row.names=row.names,fix.empty.names=F)
              print(mat,quote=FALSE)
              #
              # if there are any state equations
              if (n>0) {
                cat("\nSystem Equations:\n\n")
                lapply(private$sys.eqs,function(x) cat("\t",deparse1(x$form),"\n"))
              }
              if (m>0) {
                cat("\nObservation Equations:\n\n")
                for (i in 1:length(private$obs.eqs)) {
                  bool = private$obs.names[i] %in% private$obsvar.names
                  if (bool) {
                    cat("\t",deparse1(private$obs.eqs[[i]]$form),"+ e", "\t","e ~ N(0,",paste0(deparse1(private$obs.var[[i]]$rhs),")"),"\n")
                  } else {
                    cat("\t",deparse1(private$obs.eqs[[i]]$form),"+ e", "\t","e ~ N(0,?)","\n")
                  }
                }
              }
              if (q>0) {
                cat("\nAlgebraic Relations:\n\n")
                for(i in 1:length(private$alg.eqs)){
                  cat("\t",deparse1(private$alg.eqs[[i]]$form),"\n")
                }
              }
              # if the algebraic transformations have occured
              if (!is.null(private$sys.eqs.trans)) {
                cat("\nTransformed System Equations:\n\n")
                lapply(private$sys.eqs.trans,function(x) cat("\t",deparse1(x$form),"\n"))
              }
              # if the algebraic transformations have occured
              if (!is.null(private$obs.eqs.trans) | !is.null(private$obs.var.trans)) {
                cat("\nTransformed Observation Equations:\n\n")
                for (i in 1:length(private$obs.eqs)) {
                  bool = private$obs.names[i] %in% private$obsvar.names
                  if (bool) {
                    cat("\t",deparse1(private$obs.eqs.trans[[i]]$form),"+ e", "\t","e ~ N(0,",paste0(deparse1(private$obs.var.trans[[i]]$rhs),")"),"\n")
                  } else {
                    cat("\t",deparse1(private$obs.eqs[[i]]$form),"+ e", "\t","e ~ N(0,?)","\n")
                  }
                }
              }
              if (p>1) {
                cat("\nInputs:\n")
                cat("\t", paste(private$input.names[!private$input.names %in% "t"],collapse=", "))
              }
              if (par>0) {
                cat("\n\nParameters:\n")
                cat("\t", paste(private$parameter.names,collapse=", "))
              }
              if (fixedpars>0) {
                cat("\n\nFixed Parameters:\n")
                cat("\t", paste(names(private$fixed.pars),collapse=", "))
              }
              return(invisible(self))
            }
)

sdemTMB$set("public","summary",
            function(correlation=FALSE) {

              # check if model was estimated
              if (is.null(private$fit)) {
                message("Please estimate your model to get a fit summary.")
                return(invisible(NULL))
              }

              sumfit = summary(private$fit,correlation=correlation)
              return(invisible(sumfit$parameters))
            }
)


sdemTMB$set("public","plot",
            function(plot.obs=1, use.ggplot=FALSE, extended=FALSE, ggtheme=getggplot2theme()){

              # check if model was estimated
              if (is.null(private$fit)) {
                message("Please estimate your model in order to plot residuals.")
                return(invisible(NULL))
              }

              # if we have estimated
              plotlist = plot(private$fit, plot.obs=plot.obs, use.ggplot=use.ggplot,
                              extended=extended, ggtheme=ggtheme)
              return(invisible(plotlist))
            }
)
