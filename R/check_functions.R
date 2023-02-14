#######################################################
# CHECK EXISTENCE OF NAME FUNCTION
#######################################################

is_this_a_new_name = function(name, parameter_names) {
  trigger=TRUE
  if (name %in% parameter_names) {
    trigger=FALSE
  }
  return(trigger)
}

#######################################################
# CHECK EXISTENCE OF NAME FUNCTION
#######################################################

check_name = function(newvar, type, self, private) {
  # This function checks if e.g. a newly defined state name is not already
  # used by an observation, obs. variance or input. If the state name is
  # already used by a state, then we allow to continue - it will
  # then be overwritten in the add_systems method


  if (newvar %in% names(private$sys.eqs) & type!="state") {
    stop("The variable ", newvar, " is already in use (state).")
  }

  if (newvar %in% names(private$obs.eqs) & type!="obs" & type!="obsvar") {
    stop("The variable ", newvar, " is already in use (observation).")
  }

  if (newvar %in% names(private$inputs) & type!="input") {
    stop("The variable ", newvar, " is already in use (input).")
  }

  return(invisible(self))
}

#######################################################
# CHECK SYSTEM EQUATION FUNCTION
#######################################################

check_system_eqs = function(form, self, private, silent=FALSE) {

  # CHECK "FORM"
  #######################################################

  # The object quote(dx ~ x ) is of class "call"
  if(!inherits(form,"formula")){
    stop("The system equation should be a formula e.g. dx ~ ... * dt + ... * dw1 + ... * dw2")
  }

  lhs = form[[2]]
  rhs = form[[3]]
  # If  LHS has length one (single term) and class "name"
  if (!(length(lhs) == 1)) {
    stop("You have multiple terms on the left-hand side")
  }
  # Is the state name valid?
  state = stringr::str_match(deparse(lhs),"^d([a-zA-Z]*[a-zA-Z0-9]*)$")[2]
  if (is.na(state)) {
    stop("That state name is not allowed - use d followed by any number of letters, followed by any number of digits")
  }

  # dont use dt or dw in the state name
  match = stringr::str_match(deparse(lhs),"^(?!d[tw])")
  if (is.na(match)) {
    stop("The state name can't begin with dt or dw")
  }

  # CHECK VARIABLES IN FORM SDE
  #######################################################

  diff.proc = unique(unlist(stringr::str_extract_all(paste(deparse(rhs),collapse=""),
                                                     "dw([a-zA-Z0-9]*)")))
  # Check if there are any diffusion processes
  if(length(diff.proc)==0){
    stop("You are not allowed to specify processes without any diffusion dw(...).")
  }
  diff.processes = c("dt",diff.proc)
  diff.terms = lapply(diff.processes, function(x) { Deriv::Deriv(rhs, x, cache.exp=FALSE) })
  names(diff.terms) = diff.processes

  # Check for dt/dw cross terms
  valid = all(unlist(lapply(diff.terms, function(x) all(is.na(match(diff.processes, all.vars(x)))))))
  if (!valid) { stop("Illegal cross terms between dt and dw") }

  # Check if any parameters are outside scope like c in: f(.) * dt + g(.) * dw + c
  pars = unique(all.vars(rhs))
  pars.after = as.vector(c(diff.processes, unlist(lapply(diff.terms, all.vars))))
  if (any(is.na(match(pars, pars.after)))) {
    stop("Illegal terms outside of the drift and diffusion terms")
  }

  # PUT THIS SOMEWHERE ELSE PUT THIS SOMEWHERE ELSEPUT THIS SOMEWHERE ELSEPUT THIS SOMEWHERE ELSE
  # CHECK FOR STATE.DEP DIFF
  #######################################################

  # State dependence in the diffusion
  # valid = all(unlist(lapply(tail(diff.terms,-1), function(x) all(is.na(match(state, all.vars(x)))))))
  # if (!valid) {
  #   # If state dep is not enabled already then ask the user
  #   if (!(private$state.dep.diff)) {
  #     warning("Beware! If you choose to specify state-dependent diffusion you are on your own.
  #             You might consider performing a Lamperti Transformation z = int 1/g(x) dx which
  #             satisfies the boundary of your system i.e. kills the noise at the boundary.
  #             You should verify that the drift behaves as desired when you approach the boundary as well.
  #             Good luck.")
  #     private$state.dep.diff = TRUE
  #   }
  # }

  # is there a state with that name already?
  if (state %in% names(private$sys.eqs) & !silent) {
    # warning("Overwriting state for ", state)
  }

  # grab allvars
  bool = unique(all.vars(rhs)) %in% diff.processes
  variables = unique(all.vars(rhs))[!bool]

  # test for linearity


  # result = list(list(form=form,rhs=rhs,diff=diff.terms,allvars=variables))
  result = list(list(form=form, rhs=rhs, allvars=variables, diff=diff.processes))
  names(result) = state
  return(result)
}

#######################################################
# CHECK OBSERVATION EQUATION FUNCTION
#######################################################

check_observation_eqs <- function(form, self, private, silent=FALSE) {

  if(!inherits(form,"formula")){
    stop("The observation equation should be a formula e.g. 'y ~ ...")
  }

  lhs <- form[[2]]
  rhs <- form[[3]]

  # is there one term on the lhs?
  if (!(length(lhs) == 1)) {
    stop("You can only have one term on the left-hand side")
  }

  # you cannot observe a differential process
  obs = stringr::str_match(deparse(lhs),"^(?!d[tw])[[:alnum:]]*")
  if (is.na(obs)) {
    stop("You can't observe a differential process.")
  }

  # is there an observation with that name already?
  if (deparse(lhs) %in% names(private$obs.eqs) & !silent) {
    # warning("Overwriting observation for ", deparse(lhs))
  }

  # grab allvars
  variables = unique(all.vars(rhs))

  result = list(list(form=form,rhs=rhs,allvars=variables))
  names(result) <- deparse(lhs)
  return(result)
}

#######################################################
# CHECK OBSERVATION EQUATION FUNCTION
#######################################################

check_observation_variance_eqs <- function(form, self, private, silent=FALSE) {

  if(!inherits(form,"formula")){
    stop("The observation equation should be a formula e.g. 'y ~ ...")
  }

  lhs <- form[[2]]
  rhs <- form[[3]]
  if (!(length(lhs) == 1)) {
    stop("You can only have one term on the left-hand side")
  }

  # Is there an observation with that name?
  if (!(deparse(lhs) %in% names(private$obs.eqs))) {
    stop("Please add an observation equaiton for ", deparse(lhs), " before specifying its variance")
  }

  # is there an observation variance with that name already?
  if (deparse(lhs) %in% names(private$obs.var) & !silent) {
    # warning("Overwriting observation variance for ", deparse(lhs))
  }

  # grab allvars
  variables = unique(all.vars(rhs))

  result = list(list(form=form,rhs=rhs,allvars=variables))
  names(result) <- deparse(lhs)
  return(result)
}

#######################################################
# CHECK DATA INPUT FUNCTION
#######################################################

check_inputs <- function(input, self, private) {

  # Check for correct input class
  if (!is.name(input)) {
    stop("The inputs should be of class name i.e. use $add_inputs(a)")
  }

  # Does the input name start with dt or dw?
  valid = !is.na(stringr::str_match(input,"^(?!d[tw])[[:alnum:]]*"))
  if (!valid) {
    stop("Input names are not allowed to start with dt or dw")
  }

  # Reserved input names
  valid = !(input == "t")
  if (!valid) {
    stop("'t' is already reserved for the time vector")
  }

  result = list(list(input=input))
  names(result) = deparse(input)
  return(result)
}

#######################################################
# CHECK PARAMETERS
#######################################################

check_parameter = function(form, self, private) {


  # is formula class?
  if (!(inherits(form,"formula"))) {
    stop("You must supply a formula")
  }

  # the parameter name strings must start with a character
  parname = deparse1(form[[2]])
  bool = stringr::str_detect(parname,"^[[:alpha:]][[:alnum:]_-]*$")
  if(!bool){
    stop("The parameter name ",parname, " is not valid. The name must begin with a letter, and can only contain numerals, letters, underscore (_) and dash (-).")
  }

  # parameter name can't begin with dw or dt
  bool = stringr::str_detect(parname,"^(?!d[tw])[[:alnum:]]*")
  if(!bool){
    stop("The parameter names are not allowed to start with dt or dw")
  }



}

check_parameter_matrix <- function(parmat, self, private) {

  # is matrix class?
  if (!(inherits(parmat,"matrix"))) {
    stop("You must supply a parameter matrix.")
  }

  # has 3 columns?
  if (!ncol(parmat)==3 & !nrow(parmat)==0 ) {
    stop("The parameter matrix must have 3 columns corresponding to initial value, lower boundary and upper boundary,
         and at least 1 row")
  }

  # are parameter names supplied?
  parnames = rownames(parmat)
  if (is.null(parnames)) {
    stop("You have not supplied any parameter names. Use rownames")
  }

  # the parameter name strings must start with a character
  bool = stringr::str_detect(parnames,"^[[:alpha:]][[:alnum:]_-]*$")
  if (sum(bool) != length(bool)) {
    stop("The parameter names ",paste(parnames[!bool],collapse = ", "), " are not valid")
  }

  # parameter name can't begin with dw or dt
  bool = stringr::str_detect(parnames,"^(?!d[tw])[[:alnum:]]*")
  if (sum(bool) != length(bool)) {
    stop("Parameter names are not allowed to start with dt or dw")
  }

  result = list(parnames)
  return(result)
}

#######################################################
# CHECK ALGEBRAIC RELATIONS
#######################################################

check_algebraics = function(form, self, private) {

  if(!inherits(form,"formula")){
    stop("The algebraic relation should be a formula e.g. 'theta ~ exp(log_theta) or x ~ logit(z)")
  }
  #
  lhs = form[[2]]
  rhs = form[[3]]
  # Only single terms on LHS
  if (!(length(lhs) == 1)) {
    stop("You have multiple terms on the left-hand side")
  }

  # no algebraics for differentials (dt/dw etc)
  deparse_lhs = paste(deparse(lhs),collapse="")
  deparse_rhs = paste(deparse(rhs),collapse="")
  bool = stringr::str_match(deparse_lhs,"^(?!d[tw])[[:alnum:]]*")
  if (is.na(bool)) {
    stop("You are not allowed to redefine differential processes.")
  }
  bool = stringr::str_match(deparse_rhs,"^(?!d[tw])[[:alnum:]]*")
  if (is.na(bool)) {
    stop("No differential processes on the right-hand side.")
  }

  # Is there already an algebraic relation for this variable?
  if (deparse(lhs) %in% names(private$alg.eqs)) {
    # warning("Overwriting algebraic for ", deparse(lhs))
  }

  result = list(list(form=form,rhs=rhs))
  names(result) = deparse_lhs
  return(result)
}

#######################################################
# CHECK CONSTANTS
#######################################################

check_constants = function(form, self, private) {

  if(!inherits(form,"formula")){
    stop("The algebraic relation should be a formula e.g. 'theta ~ 2")
  }
  #
  lhs = form[[2]]
  rhs = form[[3]]
  # Only single terms on LHS
  if (!(length(lhs) == 1)) {
    stop("You have multiple terms on the left-hand side")
  }

  if (!is.numeric(rhs)) {
    stop("The right-hand side must be a numeric value")
  }

  # Is there already a constant for this variable?
  if (deparse(lhs) %in% names(private$alg.eqs)) {
    # warning("Overwriting constant ", deparse(lhs))
  }

  result = list(list(form=form,rhs=rhs))
  names(result) = deparse(lhs)
  return(result)

}

#######################################################
# CHECK MODEL CHANGES
#######################################################

was_model_already_built = function(self, private) {

  if (private$build) {
    message("You changed the model. I will recompile on next call to $estimate. Use $set_compile(FALSE) to avoid.")
    private$compile = TRUE
    private$build = FALSE
  }

  return(invisible(self))
}
