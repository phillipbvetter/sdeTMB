write_TMB_cpp = function(model,data){
  
  sdeEq = model$sdeEq
  obsEq = model$obsEq
  
  n = length(sdeEq)
  m = length(obsEq)
  # Get state and observation variables
  state = c()
  rhs = list()
  for(i in 1:n){
    state[i] = deparse(as.name(sub("^d?([[:alnum:]]+)", "\\1", sdeEq[[i]][[1]][[2]])))
    rhs[[i]]   = sdeEq[[i]][[1]][[3]]
  }
  obs = c()
  for(i in 1:length(obsEq)){
    obs[i]   = all.vars(obsEq[[i]][[1]][[2]])
  }
  
  # Get the drift and diffusion terms (dt, dw1, dw2...)
  diff.processes = c("dt",sprintf(rep("dw%i",n),1:n))
  diff.terms = list()
  for(i in 1:n){
    diff.terms[[i]]        = lapply(diff.processes, function(x) { D(rhs[[i]], x) })
    names(diff.terms[[i]]) = diff.processes
  }
  
  # Find time-dep inputs
  timedepInput = sort(unique(unlist(lapply(sdeEq,all.vars))))
  timedepInput = timedepInput[!timedepInput %in% diff.processes]
  timedepInput = timedepInput[!timedepInput %in% paste("d",state,sep="")]
  timedepInput = timedepInput[!timedepInput %in% names(c(data$pars, data$constants))]
  
  #Initialize C++ file
  full_modelname = paste(model$modelname,".cpp",sep="")
  fileconn = file(full_modelname)
  
  txt = "#include <TMB.hpp>"
  txt = append(txt, "using namespace density;")
  writeLines(txt,full_modelname)
  
  # Construct diffusion matrix function
  gvars0 = c()
  for(i in 1:n){
    gvars0 = as.vector(c(gvars0, unique(unlist(lapply(diff.terms[[i]][-1],all.vars)))))
  }
  gvars1 = paste("Type", sort(gvars0), collapse=", ")
  gvars2 = gvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      gvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=gvars2)
    }
  }
  
  temptxt = sprintf("template<class Type>\nmatrix<Type> Gmat(%s){",gvars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> G(%i,%i);",n,n))
  for(i in 1:n){
    for(j in 1:n){
      term = paste(deparse(hat2pow(diff.terms[[i]][[j+1]])),collapse="")
      temptxt = append( temptxt , sprintf("\tG(%i,%i) = %s;",i-1,j-1,term))
    }
  }
  temptxt = append(temptxt, sprintf("\tmatrix<Type> I(%s,%s);\n\tI.setIdentity();\n\tI *= 1e-4;\n\tmatrix<Type> K = G+I;",n,n))
  temptxt = append(temptxt, "\treturn(K);\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  # Construct drift function
  fvars0 = sort(unique(unlist(lapply(lapply(rhs,function(x) D(x,"dt")),all.vars))))
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars2 = fvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      fvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=fvars2)
    }
  }
  temptxt = sprintf("template<class Type>\nvector<Type> fvec(%s){",fvars1)
  temptxt = append(temptxt, sprintf("\tvector<Type> f(%i);",n))
  for(i in 1:n){
    term = paste(deparse(hat2pow(diff.terms[[i]]$dt)),collapse="")
    temptxt = append(temptxt, sprintf("\tf(%i) = %s;",i-1,term))
  }
  temptxt = append(temptxt, "\treturn(f);\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  txt = append(txt,"template<class Type>\nType objective_function<Type>::operator() ()
{")
  writeLines(txt,full_modelname)
  
  NOTdataNames = c("constants","pars",paste("obsvar",obs,sep=""),state)
  dataNames = names(data)[!(names(data) %in% NOTdataNames)]
  # Data and Parameters
  for(i in 1:length(dataNames)){
    nam = dataNames[i]
    txt = append(txt, sprintf("\tDATA_VECTOR(%s);",nam),length(txt))
  }
  writeLines(txt,full_modelname)
  # Hidden states
  for(i in 1:length(state)){
    txt = append(txt,sprintf("\tPARAMETER_VECTOR(%s);",state[i]),length(txt))
  }
  writeLines(txt,full_modelname)
  # Parameters
  for(i in 1:length(data$pars)){
    nam = names(data$pars)[i]
    txt = append(txt, sprintf("\tPARAMETER(%s);",nam))
    
  }
  writeLines(txt,full_modelname)
  # Constants
  # for(i in 1:length(data$constants)){
  #   nam = names(data$constants)[i]
  #   if(length(data$constants[[i]])>1){
  #     txt = append( txt , paste(sprintf("\tvector<Type> %s <<",nam),paste(data$constants[[i]],collapse=",")))
  #   } else {
  #     txt = append( txt , sprintf("\tType %s = %f;",nam,data$constants[[i]]))
  #   }
  # }
  for(i in 1:length(data$constants)){
    nam = names(data$constants)[i]
    if(length(data$constants[[i]])>1){
      txt = append( txt , sprintf("\tDATA_VECTOR(%s);",nam))
    } else {
      txt = append( txt , sprintf("\tDATA_SCALAR(%s);",nam))
    }
  }
  writeLines(txt,full_modelname)
  
  # Likelihood
  txt = append(txt,"\tType nll = 0;")
  writeLines(txt,full_modelname)
  
  # Create variables
  temptxt = sprintf("\tvector<Type> f(%i);",n)
  temptxt = append(temptxt, sprintf("\tvector<Type> Xi(%i);",n))
  temptxt = append(temptxt, sprintf("\tvector<Type> Xip1(%i);",n))
  temptxt = append(temptxt, sprintf("\tvector<Type> Z(%i);",n))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> G(%i,%i);",n,n))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> V(%i,%i);",n,n))
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  # dt for loop
  txt = append(txt,"\tType dt;\n\tfor(int i=0;i<t.size()-1;i++){\n\t\tdt = t(i+1) - t(i);" )
  writeLines(txt,full_modelname)
  
  txt = append(txt, sprintf("\t\tf = fvec(%s);",paste(fvars2,collapse=", ")))
  txt = append(txt, sprintf("\t\tG = Gmat(%s);",paste(gvars2,collapse=", ")))
  txt = append(txt, sprintf("\t\tXi << %s;",paste(state,"(i)",collapse=", ",sep="")))
  txt = append(txt, sprintf("\t\tXip1 << %s;",paste(state,"(i+1)",collapse=", ",sep="")))
  txt = append(txt, sprintf("\t\tZ = Xip1 - Xi - f * dt;",paste(state,"(i)",collapse=", ",sep="")))
  txt = append(txt, sprintf("\t\tV = G*G.transpose()*dt;"))
  txt = append(txt, "\t\tnll += MVNORM(V)(Z);\n\t}")
  writeLines(txt,full_modelname)
  
  #In the non-linear TMB case we need to alter the likelihood contribution from the variance because they can be coupled
  
  # Write observation contribution to likelihood
  obsvars0 = unlist(lapply(obsEq,function(x) deparse(x[[1]][[2]])))
  obsvars2 = paste(obsvars0,"(i)",sep="")
  
  hvars0 = unlist(lapply(lapply(obsEq,function(x) x[[1]][[3]]), all.vars))
  hvars0 = sort(unique(hvars0))
  hvars1 = paste("Type", hvars0, collapse=", ")
  hvars2 = hvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      hvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(j)",timedepInput[i]), x=hvars2)
    }
  }
  
  # Write observation contribution to likelihood
  for(i in 1:m){
    txt = append(txt, sprintf("\tfor(int i=0;i<%s.size(); ++i){\n\t\tint j=CppAD::Integer(%s(i));",obs[i],sprintf(rep("iobs%s",m),obs)[i]))
    txt = append(txt, sprintf("\t\tnll -= dnorm(%s,%s,%s,1);\n\t}",obsvars2[i],hvars2[i],deparse(model$obsVar[[i]][[1]])))
  }
  
  txt = append(txt, "return(nll);\n}")
  writeLines(txt,full_modelname)
  
  # Close file connection
  close(fileconn)
}
