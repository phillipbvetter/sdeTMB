checkModelAndData = function(model,data,control,method,state.dependent.diffusion){
  
  model_info = getModelInfo(model,data)
  names_model_info = names(model_info)
  for(i in 1:length(model_info)){
    assign(names_model_info[i],model_info[[i]])
  }
  
  # Check if chosen method exists
  available_methods = c("tmb","tmb_exact","ekf","ukf","kalman_adaptive")
  if(!any(is.element( method, available_methods ))){
    stop(sprintf("The chosen method is not valid. Choose one of: %s",paste(available_methods,collapse=", ")))
  }
  
  #Check if parsed model has mandatory entries
  modelnames = names(model)
  if(!any(is.element( modelnames, "modelname" ))){
    stop("Your model list does not contain an entry 'modelname'")
  }
  if(!any(is.element( modelnames, "sdeEq" ))){
    stop("Your model list does not contain an entry 'sdeEq'")
  }
  if(!any(is.element( modelnames, "obsEq" ))){
    stop("Your model list does not contain an entry 'obsEq'")
  }
  if(!any(is.element( modelnames, "obsVar" ))){
    stop("Your model list does not contain an entry 'obsVar'")
  }
  
  #Check if the time vector is increasing as expected
  if(any(diff(data$t)<=0)){
    ids = which(diff(data$t)<=0)
    stop(sprintf("The time-vector is non-increasing at the following indice(s) %s",paste(ids,collapse=", ")))
  }
  
  # State-dependent diffusion is not allowed, only in "playground" mode
  if(!state.dependent.diffusion){
    n.diff = length(diff.processes)-1
    for(i in 1:n){
      for(j in 1:n.diff){
        eq = diff.terms[[i]][[j+1]]
        prev.vars = all.vars(eq)
        present.states = state[state %in% prev.vars]
        if( length(present.states) > 0 ){
          stop("The system has state-dependent diffusion, which is not well-supported. 
               Consider performing a Lamperti transformation to achieve state-independence. 
               If you wish to continue regardless, parse 'state.dependent.diffusion=TRUE' to MakeNLL")
        }
      }
    } 
  }
  
  # Find time-dep inputs for later checks
  timedepInput = sort(unique(unlist(lapply(sdeEq,all.vars))))
  timedepInput = timedepInput[!timedepInput %in% diff.processes]
  timedepInput = timedepInput[!timedepInput %in% paste("d",state,sep="")]
  timedepInput = timedepInput[!timedepInput %in% names(c(data$pars, data$constants))]
  
  # Check data-list has mandatory entries for the different methods
  data.names = names(data)
  switch(method,
         tmb = {
           if(!any(is.element( data.names, "t" ))){
             stop("Your data-list does not contain a time vector t")
           }
           if(!all(is.element(state,data.names))){
             str = state[!is.element(state,data.names)]
             stop(sprintf("Your data-list does not contain initial values for the random effects vector(s): %s",paste(str,collapse=", ")))
           }
           if(!any(is.element( data.names, "pars" ))){
             stop("Your data-list does not contain a 'pars' entry with initial parameter values")
           }
           if(!all(sapply(obs,function(x) any(is.element(data.names,x))))){
             bool = !sapply(obs,function(x) any(is.element(data.names,x)))
             str = paste(names(bool)[bool],collapse=", ")
             stop(sprintf("The following observed variables are not present in the data-list: %s", str ))
           }
           if(length(timedepInput)>0){
             for(i in 1:length(timedepInput)){
               if(any(is.na(data[[timedepInput[i]]]))){
                 ids = which(is.na(data[[timedepInput[i]]]))
                 stop(sprintf("The input vector %s has NA-values at the indice(s) %s",timedepInput[i],paste(ids,collapse=", ")))
               }
             }
           }
         },
         ekf =,
         ukf = {
           # does time vector exist?
           if(!any(is.element( data.names, "t" ))){
             stop("Your data-list does not contain a time vector t")
           }
           # does initial state vector X0 exist?
           if(!any(is.element( data.names, "X0" ))){
             stop("Your data-list does not contain an initial state vector X0")
           }
           # does X0 have correct dimensions?
           if(!length(data[["X0"]])==n){
             stop(sprintf("Your initial state X0 does not have length %s",n))
           }
           # does initial state-covariance matrix P0 exist?
           if(!any(is.element( data.names, "P0" ))){
             stop("Your data-list does not contain an initial state covariance matrix P0")
           }
           # does P0 have correct dimensions?
           if(!all(dim(data[["P0"]]) == c(n,n))){
             stop(sprintf("Your initial state covariance P0 does not have dimensions [%s,%s]",n,n))
           }
           # does time-step scalar dt exist?
           if(!any(is.element( data.names, "dt" ))){
             stop("Your data-list does not contain a time-step dt for the forward-euler scheme of the moment differential equations")
           }
           # does parameter list pars exist?
           if(!any(is.element( data.names, "pars" ))){
             stop("Your data-list does not contain a 'pars' entry with initial parameter values")
           }
           # does data for all observed variables exist?
           if(!all(sapply(obs,function(x) any(is.element(data.names,x))))){
             bool = !sapply(obs,function(x) any(is.element(data.names,x)))
             str = paste(names(bool)[bool],collapse=", ")
             stop(sprintf("The following observed variable names are not present in the data-list: %s", str ))
           }
           # does input vectors have any NA values?
           timedepInput = timedepInput[!timedepInput %in% state]
           if(length(timedepInput)>0){
             for(i in 1:length(timedepInput)){
               if(any(is.na(data[[timedepInput[i]]]))){
                 ids = which(is.na(data[[timedepInput[i]]]))
                 stop(sprintf("The input vector %s has NA-values at the indice(s) %s",timedepInput[i],paste(ids,collapse=", ")))
               }
             }
           }
         }
  )
  
  # Check the MAP mean and covariance matrices
  if(is.null(data$MAPbool)){
    data$MAPbool = 0
  }
  is.scalar.zeroORone <- function (x) length(x) == 1L && is.vector(x, mode = "numeric") && (x==0 || x==1)
  if(!is.scalar.zeroORone(data$MAPbool)){
    stop("The data$MAPbool value can only be 0 or 1. Select 1 to enable MAP estimation.")
  }
  if(data$MAPbool==1){
    if(!is.vector(data$MAPmean,mode="numeric")){
      stop("Please provide a valid prior mean vector for MAP estimation.")
    }
    if(!is.vector(as.vector(data$MAPcov),mode="numeric")){
      stop("Please provide a valid prior covariance matrix for MAP estimation.")
    }
    if(!(length(data$MAPmean)==length(data$pars))){
      stop(sprintf("The MAP mean vector has %s entries but should have %s",length(data$MAPmean),length(data$pars)))
    }
    if(!(length(data$MAPcov) == length(data$pars)^2)){
      stop(sprintf("The MAP covariance matrix has length %s, but should have length %s",length(data$MAPcov),length(data$pars)^2))
    }
    if(is.null(dim(data$MAPcov))){
      stop("The provided MAP covariance matrix is not a matrix, but a vector.")
    }
    if(!(all(dim(data$MAPcov)==rep(length(data$pars),2)))){
      stop(sprintf("The MAP covariance matrix has dimensions (%s) but should have (%s)",paste(dim(data$MAPcov),collapse=","),paste(rep(length(data$pars),2),collapse=",")))
    }
  }
  
  # Check if all variables in SDE formulation are accounted for
  allvars = c()
  for(i in 1:n){
    allvars = c(allvars,unlist(sapply(diff.terms[[i]],all.vars)))
  }
  allvars = unique(allvars)
  allvars = allvars[!(allvars %in% state)]
  allvars = allvars[!(allvars %in% data.names)]
  allvars = allvars[!(allvars %in% names(data$pars))]
  allvars = allvars[!(allvars %in% names(data$constants))]
  if(length(allvars)>0){
    stop(sprintf("The following entries are not accounted for: %s",paste(allvars,collapse=", ")))
  }
  
  # check control list
  if(!is.null(control$loss)){
    if(!(control$loss %in% c("identity","tukey","huber"))){
      stop(sprintf("The chosen loss regularization %s does not exist. Choose either identity, tukey or huber"),control$loss)
    }
  }
  is.positive.scalar <- function (x) length(x) == 1L && is.vector(x, mode = "numeric") && x>0
  if(!is.positive.scalar(control$c)){
    stop(sprintf("Please choose a positive integer value for control$c"))
  }
  
  return(data)
}

