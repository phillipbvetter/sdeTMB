#' @title Constructor for R6 ctsmrTMB class
#' @description Call ctsmrTMB$new() to construct a new instance of the class object
#' for specifying a new stochastic state space model.
#' @returns A model object of class 'ctsmrTMB'
#' @export
ctsmrTMB = R6::R6Class(
  classname = "ctsmrTMB",
  # public fields
  public = list(
    #' @description Construct new ctsmrTMB model object
    initialize = function() {
      message("Created New ctsmrTMB Object")
      # public fields
      #
      # modelname
      private$modelname = "sde_model"
      private$cppfile.directory = getwd()
      private$cppfile.path = paste(getwd(),"/",private$modelname,sep="")
      private$cppfile.path.with.method = NULL
      private$modelname.with.method = NULL
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
      private$control.nlminb = list()
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
    ########################################################################
    ########################################################################
    #' @description 
    #' Define and add multiple stochastic differential equation governing the process of individual state variables 
    #' on the form 
    #' 
    #' \code{d<state> ~ f(t,<states>,<inputs>) * dt + g1(t,<states>,<inputs>) * dw1 
    #'                                          + g2(t,<states>,<inputs>) * dw2 
    #'                                            + ...}
    #'                                            
    #' where \code{f} is the drift, and \code{g1, g2, ...} are diffusions, with differential brownian motions dw1, dw2, ...
    #' 
    #' @examples 
    #' # Specify Ornstein-Uhlenbeck Process
    #' add_systems(dX ~ theta * (mu - X + u) * dt + sigma * dw)
    #' 
    #' # Specify Lokta-Volterra System of Equations
    #' add_systems( dN ~ ( r * N * (1-N/K) - c*N*P/(N+Nbar) ) * dt + sigmaN * N * dw1,
    #'              dP ~ ( eps*c*N*P/(N+Nbar) - mu * P ) * dt + sigmaP * P * dw2 )
    #' @param form formula specifying the stochastic differential equation to be added to the system.
    #' @param ... additional formulas similar to \code{form} for specifying multiple equations at once.
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
      # if model was already built, set compile flag to true
      was_model_already_built(self, private)
      #
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    #' @description
    #' Define and add a relationship between an observed variable and system states. The observation equation
    #' takes the form
    #' 
    #' \code{<observation> ~ h(t,<states>,<inputs>) + e)}
    #' 
    #' where \code{h} is the observation function, and \code{e} is normally distributed noise with zero mean and variance
    #' to be specified. The observation variable should be present in the data provided when calling
    #' \code{estimate(.data)} for parameter estimation.
    #'  
    #' @examples
    #' #Specify observation directly as a latent state
    #' add_observations(y ~ x)
    #' 
    #' Specify observation as the sum of exponentials of two latent states
    #' add_observations(y ~ exp(x1) + exp(x2))
    #' @param form formula class specifying the obsevation equation to be added to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple observation equations at a time.
    add_observations = function(form,...) {
      # lapply for multi-arguments
      lapply(c(form,...), function(form) {
        # Function Body
        res = check_observation_eqs(form, self, private)
        check_name(names(res), "obs", self, private)
        private$obs.eqs[[names(res)]] = res[[1]]
        private$obs.names = names(private$obs.eqs)
        # Function Body
      })
      # if model was already built, set compile flag to true
      was_model_already_built(self, private)
      #
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    #' @description Specify the variance of an observation equation.
    #' 
    #' A defined observation variable \code{y} in e.g. \code{add_observations(y ~ 
    #' h(t,<states>,<inputs>)} is pertubed by Gaussian noise with zero mean and variance 
    #' to-be specified using \code{add_observation_variances(y ~ p(t,<states>,<inputs>)}. We can
    #' for instance declare \code{add_observation_variances(y ~ sigma_x^2} where \code{sigma_x} 
    #' is a fixed effect parameter to be declared through \code{add_parameters}.
    #' 
    #' @param form formula class specifying the obsevation equation to be added to the system.
    #' @param ... additional formulas identical to \code{form} to specify multiple observation equations at a time.
    add_observation_variances = function(form,...) {
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
    },
    ########################################################################
    ########################################################################
    #' @description Declare variables as data inputs
    #'
    #' Declare whether a variable contained in system, observation or observation variance equations 
    #' is an input variable. If e.g. the system equation contains an input variable \code{u} then it
    #' is declared using \code{add_inputs(u)}. The input data \code{u} must be contained in the 
    #' data.frame \code{data} when calling \code{estimate(data)} for parameter estimation.
    #' 
    #' @param ... additional formulas identical to \code{form} to specify multiple observation equations at a time.
    add_inputs =  function(...) {
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
    },
    ########################################################################
    ########################################################################
    #' @description Declare variables as fixed effects and specify their initial value, lower and
    #' upper bound used when calling the maximum likelihood optimization.
    #' 
    #' There are two ways to declare parameters. You can declare parameters using formulas i.e.
    #' \code{add_parameters( theta ~ c(1,0,10), mu ~ c(0,-10,10) )}, where the values are initial,
    #' lower and upper bound respectively. Alternatively you can provide a 3-column matrix where
    #' rows corresponds to different parameters, and the parameter names are provided as rownames
    #' of the matrix.
    #'
    #' @param ... formula whose left-hand side is the parameter name, and right hand side is a vector
    #' of length 3 with inital value, lower and upper bound respectively. You can provide multiple
    #' parameters at once by seperating formulas with comma.
    #' @param parameter.matrix matrix of 3 columns where rows correspond to variables. The variable
    #' names must be provided as rownames to the matrix. The columns are initial value, lower and 
    #' upper bound respectively.
    add_parameters = function(...,parameter.matrix=NULL) {
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
    },
    ########################################################################
    ########################################################################
    #' @description Add algebraic relations.
    #' 
    #' Algebraic relations is a convenient way to transform parameters in your equations.
    #' In the Ornstein-Uhlenbeck process the rate parameter \code{theta} is always positive, so
    #' estimation in the log-domain is a good idea. Instead of writing \code{exp(theta)} directly
    #' in the system equation one can transform into the log domain using the algebraic relation
    #' \code{add_algebraics(theta ~ exp(logtheta))}. All instances of \code{theta} is replaced
    #' by \code{exp(logtheta)} when compiling the C++ function. Note that you must provide values
    #' for \code{logtheta} now instead of \code{theta} when declaring parameters through 
    #' \code{add_parameters}
    #' 
    #' @param form formula specifying the stochastic differential equation(s) to be added to the system.
    #' @param ... additional formulas similar to \code{form} for specifying multiple equations at once.
    add_algebraics = function(form,...) {
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
    },
    ########################################################################
    ########################################################################
    #' @description Declare variables as scalar constants.
    #'
    #' @param form formula whose left-hand side is the constant variable name, and the right-hand side
    #' is its numeric value.
    #' @param ... additional formulas similar to \code{form} for specifying multiple constants at once.
    add_constants = function(form,...) {
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
    },
    ########################################################################
    ########################################################################
    #' @description Declare the initial state values i.e. mean and covariance for the system states.
    #' 
    #' @param mean numeric vector of size equal to \code{n} (number of system states)
    #' @param cov matrix (symmetric positive semi-definite) of dimensions \code{n^2}. 
    set_initial_state = function(mean,cov) {
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
      # convert scalar to matrix
      if(!is.matrix(cov) & is.numeric(cov) & length(cov)==1){
        cov = cov*diag(1)
      }
      if (!is.numeric(cov)) {
        stop("The covariance matrix is not a numeric")
      }
      if (any(is.na(cov))) {
        stop("The covariance matrix contains NAs")
      }
      if (any(eigen(cov)$values < 0)){
        stop("The covariance matrix is not positive semi-definite")
      }
      if (!isSymmetric.matrix(cov)){
        stop("The covariance matrix is symmetric")
      }
      private$initial.state = list(mean=mean,cov=as.matrix(cov))
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    #' @description Set a Lamperti Transformation
    #'
    #' If the provided system equations have state dependent diffusion in of a few available ways
    #' then it is advantageous to perform a transformation to remove the state dependence. This 
    #' comes at the cost of a more complicated drift function. The following types of state-dependence
    #' is currently supported
    #'
    #' 0. If no transformation is desired choose 'identity' (default).
    #'
    #' 1. If the diffusion term is proportional to x * dw, then a log-transform is available
    #' 
    #' 2. If the diffusion term is proportional to x * (1-x) * dw, then a logit-transform is available
    #' 
    #' 3. If the diffusion term is proportional to sqrt(x * (1-x)) * dw, then a sqrt-logit-transform 
    #' is available
    #' 
    #' 4. If the diffusion term is proportional to x * (1-x^a) * dw, for a>0 then a power-logit
    #' transform is available
    #' 
    #' @param transform character vector - one of either "identity, "log", "logit", "sqrt-logit", or "power-logit"
    #' @param states a vector of the state names for which the chosen transformation is desired. The
    #' default (NULL) is to apply the transformation to all state equations.
    set_lamperti = function(transform,states=NULL) {
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
    },
    ########################################################################
    ########################################################################
    #' @description Set modelname used to create the C++ file for TMB
    #'
    #' When calling \code{TMB::MakeADFun} the (negative log) likelihood function is created in the
    #' directory specified by the \code{set_cppfile_directory} method with name \code{<modelname>.cpp}
    #' 
    #' @param name string defining the model name.
    set_modelname = function(name) {
      # was a string passed?
      if (!is.character(name)) {
        stop("The modelname must be a string")
      }
      private$modelname = name
      private$cppfile.path = paste(private$cppfile.directory,"/",name,sep="")
    },
    ########################################################################
    ########################################################################
    #' @description Set the path directory where the constructed C++ file is created.
    #'
    #' @param directory string specifying the local path / directory
    set_cppfile_directory = function(directory) {
      # was a string passed?
      if (!is.character(directory)) {
        stop("You must pass a string")
      }
      # create directory if it does not exist?
      if (!dir.exists(directory)) {
        message("The specified directory does not exist - creating it")
        dir.create(directory)
      }
      private$cppfile.directory = directory
      
      # update private$cppfile.path by calling set_modelname
      self$set_modelname(private$modelname)
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    #' @description Enable maximum a posterior (MAP) estimation.
    #'
    #' Adds a maximum a posterior contribution to the (negative log) likelihood 
    #' function by evaluating the fixed effects parameters in a multivariate Gaussian 
    #' with \code{mean} and \code{covariance} as provided.
    #' 
    #' @param mean mean vector of the Gaussian prior parameter distribution
    #' @param cov covariance matrix of the Gaussian prior parameter distribution
    set_map = function(mean,cov) {
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
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    #' @description Construct and extract function handlers for the negative
    #' log likelihood function.
    #'
    #' The handlers from \code{TMB}'s \code{MakeADFun} are constructed and returned. 
    #' This enables the user to e.g. determine their own optimization algorithm to use
    #' for estimation.
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.  
    #' @param ode.timestep numeric value. Sets the time step-size in numerical filtering schemes. 
    #' The defined step-size is used to calculate the number of steps between observation time-points as 
    #' defined by the provided \code{data}. If the calculated number of steps is larger than N.01 where N 
    #' is an integer, then the time-step is reduced such that exactly N+1 steps is taken between observations  
    #' The step-size is used in the two following ways depending on the
    #' chosen method:
    #' 1. Kalman filters: The time-step is used as the step-size in the
    #' numerical Forward-Euler scheme to compute the prior state mean and
    #' covariance estimate as the final time solution to the first and second
    #' order moment differential equations.
    #' 2. TMB method: The time-step is used as the step-size in the Euler-Maruyama
    #' scheme for simulating a sample path of the stochastic differential equation,
    #' which serves to link together the latent (random effects) states.
    #' @param silence boolean value. Sets the tracing information for \code{TMB::MakeADFun} in the
    #' argument \code{silent} which disables outputs from the optimization algoritm during runtime.
    #' @param compile boolean value. The default (\code{FALSE}) is to not compile the C++ objective
    #' function but assume it is already compiled and corresponds to the specified model object. It is
    #' the user's responsibility to ensure correspondence between the specified model and the precompiled
    #' C++ object. If a precompiled C++ object is not found in the specified directory i.e. 
    #' in \code{<cppfile_directory>/<modelname>/(dll/so)} then the compile flag is set to \code{TRUE}.
    #' If the user makes changes to system equations, observation equations, observation variances, 
    #' algebraic relations or lamperi transformations then the C++ object should be recompiled.
    #' @param method character vector - one of either "ekf", "ukf" or "tmb". Sets the estimation 
    #' method. The package has three available methods implemented:
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residauls as is natural 
    #' when the Gaussian (negative log) likelihood is evaluated, but if the tails of the 
    #' distribution is considered too small i.e. outliers are weighted too much, then one 
    #' can choose loss functions that accounts for this. The three available types available:
    #' 
    #' 1. Quadratic loss (\code{quadratic}).
    #' 2. Quadratic-Linear (\code{huber})
    #' 3. Quadratic-Constant (\code{tukey})
    #' 
    #' The cutoff for the Huber and Tukey loss functions are determined from a provided cutoff 
    #' parameter \code{loss_c}. The implementations of these losses are approximations (pseudo-huber and sigmoid 
    #' approxmation respectively) for smooth derivatives.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    construct_nll = function(data,
                       ode.timestep=NULL, 
                       silence=FALSE, 
                       compile=FALSE,
                       method="ekf",
                       loss="quadratic",
                       loss_c=3) {
      
      # set flags
      private$set_timestep(ode.timestep)
      private$set_silence(silence)
      private$set_compile(compile)
      private$set_method(method)
      private$set_loss(loss,loss_c)
      
      # build model
      message("Building model...")
      private$build_model()
      
      # check and set data
      message("Checking data...")
      check_and_set_data(data, self, private)
      
      # construct neg. log-likelihood
      message("Constructing objective function...")
      construct_nll(self, private)
      
      # return
      message("Succesfully returned function handlers")
      return(private$nll)
    },
    ########################################################################
    ########################################################################
    #' @description Estimate the fixed effects parameters in the specified model.
    #' 
    #' @param data data.frame containing time-vector 't', observations and inputs. The observations
    #' can take \code{NA}-values.  
    #' @param use.hessian boolean value. The default (\code{TRUE}) causes the optimization algorithm
    #' \code{stats::nlminb} to use the fixed effects hessian of the (negative log) likelihood when
    #' performing the optimization. This feature is only available for the kalman filter methods 
    #' without any random effects.
    #' @param ode.timestep numeric value. Sets the time step-size in numerical filtering schemes. 
    #' The defined step-size is used to calculate the number of steps between observation time-points as 
    #' defined by the provided \code{data}. If the calculated number of steps is larger than N.01 where N 
    #' is an integer, then the time-step is reduced such that exactly N+1 steps is taken between observations  
    #' The step-size is used in the two following ways depending on the
    #' chosen method:
    #' 1. Kalman filters: The time-step is used as the step-size in the
    #' numerical Forward-Euler scheme to compute the prior state mean and
    #' covariance estimate as the final time solution to the first and second
    #' order moment differential equations.
    #' 2. TMB method: The time-step is used as the step-size in the Euler-Maruyama
    #' scheme for simulating a sample path of the stochastic differential equation,
    #' which serves to link together the latent (random effects) states.
    #' @param silence boolean value. Sets the tracing information for \code{TMB::MakeADFun} in the
    #' argument \code{silent} which disables outputs from the optimization algoritm during runtime.
    #' @param compile boolean value. The default (\code{FALSE}) is to not compile the C++ objective
    #' function but assume it is already compiled and corresponds to the specified model object. It is
    #' the user's responsibility to ensure correspondence between the specified model and the precompiled
    #' C++ object. If a precompiled C++ object is not found in the specified directory i.e. 
    #' in \code{<cppfile_directory>/<modelname>/(dll/so)} then the compile flag is set to \code{TRUE}.
    #' If the user makes changes to system equations, observation equations, observation variances, 
    #' algebraic relations or lamperi transformations then the C++ object should be recompiled.
    #' @param method character vector - one of either "ekf", "ukf" or "tmb". Sets the estimation 
    #' method. The package has three available methods implemented:
    #' 1. The natural TMB-style formulation where latent states are considered random effects
    #' and are integrated out using the Laplace approximation. This method only yields the gradient
    #' of the (negative log) likelihood function with respect to the fixed effects for optimization.
    #' The method is slower although probably has some precision advantages, and allows for non-Gaussian
    #' observation noise (not yet implemented). One-step / K-step residuals are not yet available in
    #' the package.
    #' 2. (Continous-Discrete) Extended Kalman Filter where the system dynamics are linearized
    #' to handle potential non-linearities. This is computationally the fastest method.
    #' 3. (Continous-Discrete) Unscented Kalman Filter. This is a higher order non-linear Kalman Filter
    #' which improves the mean and covariance estimates when the system display high nonlinearity, and
    #' circumvents the necessity to compute the jacobian of the drift and observation functions.
    #' 
    #' All package features are currently available for the kalman filters, while TMB is limited to
    #' parameter estimation. In particular, it is straight-forward to obtain k-step-ahead predictions
    #' with these methods (use the \code{predict} S3 method), and stochastic simulation is also available 
    #' in the cases where long prediction horizons are sought, where the normality assumption will be 
    #' inaccurate.
    #' @param loss character vector. Sets the loss function type (only implemented for the kalman filter
    #' methods). The loss function is per default quadratic in the one-step residauls as is natural 
    #' when the Gaussian (negative log) likelihood is evaluated, but if the tails of the 
    #' distribution is considered too small i.e. outliers are weighted too much, then one 
    #' can choose loss functions that accounts for this. The three available types available:
    #' 
    #' 1. Quadratic loss (\code{quadratic}).
    #' 2. Quadratic-Linear (\code{huber})
    #' 3. Quadratic-Constant (\code{tukey})
    #' 
    #' The cutoff for the Huber and Tukey loss functions are determined from a provided cutoff 
    #' parameter \code{loss_c}. The implementations of these losses are approximations (pseudo-huber and sigmoid 
    #' approxmation respectively) for smooth derivatives.
    #' @param loss_c cutoff value for huber and tukey loss functions. Defaults to \code{c=3}
    #' @param control list of control parameters parsed to \code{nlminb} as its \code{control} argument. 
    #' See \code{?stats::nlminb} for more information
    estimate = function(data, 
                        use.hessian=FALSE,
                        ode.timestep=NULL, 
                        silence=TRUE, 
                        compile=FALSE,
                        method="ekf",
                        loss="quadratic",
                        loss_c=3,
                        control=list(trace=1)) {
      
      # set flags
      private$use_hessian(use.hessian)
      private$set_timestep(ode.timestep)
      private$set_silence(silence)
      private$set_compile(compile)
      private$set_method(method)
      private$set_loss(loss,loss_c)
      private$set_control(control)
      
      
      # build model
      # if (private$compile) {
      message("Building model...")
      private$build_model()
      # }
      
      # check and set data
      message("Checking data...")
      check_and_set_data(data, self, private)
      
      # construct neg. log-likelihood function
      message("Constructing objective function...")
      construct_nll(self, private)
      
      # estimate
      message("Minimizing the negative log-likelihood...")
      optimise_nll(self, private)
      
      # create return fit
      create_return_fit(self, private)
      
      # return fit
      return(invisible(private$fit))
    },
    ########################################################################
    ########################################################################
    # PRINT FUNCTION
    #' @description Function to print the model object
    print = function() {
      n = length(private$sys.eqs)
      m = length(private$obs.eqs)
      p = length(private$inputs)-1
      ng = max(length(unique(unlist(lapply(private$sys.eqs, function(x) x$diff))))-1,0)
      q = length(private$alg.eqs)
      par = length(private$parameters)
      fixedpars = length(private$fixed.pars)
      # If the model is empty
      cat("Stochastic State Space Model:")
      basic.data = c(private$modelname,n,ng,m,p,par)
      row.names = c("Name", "States","Diffusions",
                    "Observations","Inputs",
                    "Parameters")
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
    },
    ########################################################################
    ########################################################################
    # SUMMARY FUNCTION
    #' @description Summary function for fit
    #' @param correlation boolean value. The default (\code{FALSE}) is to not provide the fixed effects parameter
    #' correlation matrix. 
    summary = function(correlation=FALSE) {
      
      # check if model was estimated
      if (is.null(private$fit)) {
        message("Please estimate your model to get a fit summary.")
        return(invisible(NULL))
      }
      
      sumfit = summary(private$fit,correlation=correlation)
      return(invisible(sumfit$parameters))
    },
    ########################################################################
    ########################################################################
    # PLOT FUNCTION
    #' @description Function to print the model object
    #' @param plot.obs integer value. The value determines which state residual plot is shown, when the states are numbered
    #' by the order in which they were defined in the model object.
    #' @param use.ggplot boolean value. The default \code{FALSE} is to use base plots in R, if \code{TRUE} then \code{ggplot2} 
    #' is used to show and return a list of all residual state plots.
    #' @param extended boolean value. No functionality yet.
    #' @param ggtheme ggplot theme. This determines the theme used to construct the ggplot2 if \code{use.ggplot=TRUE}.
    plot = function(plot.obs=1, use.ggplot=FALSE, extended=FALSE, ggtheme=getggplot2theme()){
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
  ),
  # private fields
  private = list(
    # modelname
    modelname = character(0),
    cppfile.directory = NULL,
    cppfile.path = character(0),
    cppfile.path.with.method = NULL,
    modelname.with.method = NULL,
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
    control.nlminb = NULL,
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
    fit = NULL,
    ########################################################################
    ########################################################################
    # lamperti transform functions
    add_trans_systems = function(form) {
      res = check_system_eqs(form, self, private, silent=T)
      private$sys.eqs.trans[[names(res)]] = res[[1]]
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    # lamperti transform functions
    add_trans_observations = function(form) {
      res = check_observation_eqs(form, self, private, silent=T)
      private$obs.eqs.trans[[names(res)]] = res[[1]]
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    # lamperti transform functions
    add_trans_observation_variances = function(form) {
      res = check_observation_variance_eqs(form, self, private, silent=T)
      private$obs.var.trans[[names(res)]] = res[[1]]
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    # build function
    build_model = function() {
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
    },
    ########################################################################
    ########################################################################
    # compile function
    set_compile = function(bool) {
      # is bool logical
      if (!is.logical(bool)) {
        stop("You must pass a logical value")
      }
      private$compile = bool
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    # silent function
    set_silence = function(bool) {
      # is bool logical
      if (!is.logical(bool)) {
        stop("You must pass a logical value")
      }
      private$silent = bool
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    # set method
    set_method = function(method) {
      # is the method a string?
      if (!(is.character(method))) {
        stop("You must pass a string")
      }
      # choose one of the available methods
      available_methods = c("ekf","ukf","tmb")
      if (!(method %in% available_methods)) {
        stop("That method is not available. Please choose one of the following instead: \n
                  1. Extended Kalman Filter - method = 'ekf'
                  2. Unscented Kalman Filter - method = 'ukf'
                  3. Laplace Approx using Random Effects - method = 'tmb'")
      }
      private$method = method
      return(invisible(self))
    },
    ########################################################################
    ########################################################################
    # set time-step
    set_timestep = function(dt) {
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
    },
    ########################################################################
    ########################################################################
    # USE HESSIAN FUNCTION
    use_hessian = function(bool) {
      if (!is.logical(bool)) {
        stop("This must be a logical")
      }
      private$use.hessian = bool
    },
    ########################################################################
    ########################################################################
    # SET LOSS FUNCTION
    set_loss = function(loss,c=3) {
      # is the method a string?
      if (!(is.character(loss))) {
        stop("You must pass a string")
      }
      # choose one of the available methods
      available_losses = c("quadratic","huber","tukey")
      if (!(loss %in% available_losses)) {
        stop("That method is not available. Please choose one of the following instead: \n
                  1. Quadratic loss = 'quadratic'
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
          # res = sum((funtukey(rtukey) - .pars[4]*(sigmoid(rtukey,a=.pars[1],b=.pars[2])+.pars[3]))^2)
          res = sum((funtukey(rtukey) - .pars[4]*(invlogit2(rtukey,a=.pars[1],b=.pars[2])+.pars[3]))^2)
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
    },
    ########################################################################
    ########################################################################
    # SET CONTROL FOR NLMINB
    set_control = function(control) {
      # is the control a list?
      if (!(is.list(control))) {
        stop("The control argument must be a list")
      }
      private$control.nlminb = control
      return(invisible(self))
    }
  )
)
