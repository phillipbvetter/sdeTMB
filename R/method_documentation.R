# This document contains empty functions with names similar to the methods in sdeTMB
# This allow users to easily find method documentation using ?add_systems, instead of having
# to manually retrieve them from ?sdeTMB and then finding the method among the large list of other methods


########################################################################
# ADD SYSTEMS
########################################################################
#' Adding system state equations
#' 
#' @description 
#' Define and add multiple stochastic differential equation governing the 
#' process of individual state variables on the form 
#' 
#' \code{d<state> ~ f(t,<states>,<inputs>) * dt + g1(t,<states>,<inputs>) * dw1 
#' + g2(t,<states>,<inputs>) * dw2 + ... + gN(t,<states>,<inputs>) * dwN}
#'                                            
#' where \code{f} is the drift, and \code{g1, g2, ..., gN} are diffusions, with 
#' differential brownian motions dw1, dw2, ..., dwN.
#' 
#' @examples 
#' # Specify Ornstein-Uhlenbeck Process
#' add_systems(dx ~ theta * (mu - x + u) * dt + sigma * dw)
#'              
#' @param form formula specifying the stochastic differential equation to be 
#' added to the system.
#' @param ... additional formulas similar to \code{form} for specifying 
#' multiple equations at once.
#' 
add_systems = function(form,...) {
  NULL
}

########################################################################
# ADD OBSERVATIONS
########################################################################
#' Adding observation equations
#' 
#' @description
#' Define and add a relationship between an observed variable and system states. 
#' The observation equation takes the form
#' 
#' \code{<observation> ~ h(t,<states>,<inputs>) + e)}
#' 
#' where \code{h} is the observation function, and \code{e} is normally 
#' distributed noise with zero mean and variance to be specified. The 
#' observation variable should be present in the data provided when calling
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
#' @param obsnames character vector specifying the name of the observation. When the observation left-hand side
#' consists of more than just a single variable name (when its class is 'call' instead of 'name') it will be 
#' given a name on the form obs__# where # is a number, unless obsnames is provided.
add_observations = function(form,...,obsnames=NULL) {
  NULL
}

########################################################################
# ADD OBSERVATION VARIANCES
########################################################################
#' Adding observation variances
#' 
#' @description Specify the variance of an observation equation.
#' 
#' A defined observation variable \code{y} in e.g. \code{add_observations(y ~ 
#' h(t,<states>,<inputs>)} is pertubed by Gaussian noise with zero mean and 
#' variance 
#' to-be specified using \code{add_observation_variances(y ~ p(t,<states>,<inputs>)}. 
#' We can for instance declare \code{add_observation_variances(y ~ sigma_x^2} 
#' where \code{sigma_x} is a fixed effect parameter to be declared through 
#' \code{add_parameters}.
#' 
#' @param form formula class specifying the obsevation equation to be added 
#' to the system.
#' @param ... additional formulas identical to \code{form} to specify multiple 
#' observation equations at a time.
#' 
add_observation_variances = function(form,...) {
  NULL
}

########################################################################
# ADD INPUTS
########################################################################
#' Declare variables as data inputs
#' 
#' @description Declare variables as data inputs
#' 
#' Declare whether a variable contained in system, observation or observation 
#' variance equations is an input variable. If e.g. the system equation contains 
#' an input variable \code{u} then it is declared using \code{add_inputs(u)}. 
#' The input \code{u} must be contained in the data.frame \code{.data} provided 
#' when calling the \code{estimate} or \code{predict} methods.
#' 
#' @param ... variable names that specifies the name of input variables in the defined system.
#' 
add_inputs =  function(...) {
  NULL
}

########################################################################
# ADD PARAMETERS
########################################################################
#' Declare parameters and specify initial values and lower/upper bounds for 
#' optimization procedure.
#' 
#' @description Declare which variables that are (fixed effects) parameters in
#' the specified model, and specify the initial optimizer value, as well as
#' lower / upper bounds during optimization. There are two ways to declare parameters:
#' 
#' 1. You can declare parameters using formulas i.e. \code{add_parameters( 
#' theta = c(1,0,10), mu = c(0,-10,10) )}. The first value is the initial 
#' value for the optimizer, the second value is the lower optimization 
#' bound and the third value is the upper optimization bound. 
#' 
#' 2. You can provide a 3-column matrix where rows corresponds to different 
#' parameters, and the parameter names are provided as rownames of the matrix. 
#' The columns values corresponds to the description in the vector format above.
#'
#' @param ... a named vector or matrix as described above.
add_parameters = function(...) {
  NULL
}

########################################################################
# ADD ALGEBRAICS
########################################################################
#' Add algebraic relationships 
#'
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
  NULL
}

########################################################################
# SET INITIAL STATE
########################################################################
#' Set initial state mean and covariance 
#'
#' @description Declare the initial state values i.e. mean and covariance for the system states.
#' 
#' @param initial.state a named list of two entries 'x0' and 'p0' containing the initial state and covariance of the state
#' @param estimate boolean value which indicates whether or not the initial conditions
#' shall be estimated as fixed effects parameters. The provided mean and covariance are then
#' used as initial guesses
#' 
set_initial_state = function(initial.state, estimate=FALSE) {
  NULL
}

