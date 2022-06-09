write_linearExact_cpp = function(model,data){
  
  # Substitute algebraic expressions
  for(i in 1:length(model$algeqs)){
    curlist = list()
    curlist[[names(model$algeqs)[i]]] = model$algeqs[[i]][[1]]
    model$sdeEq = lapply(model$sdeEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsEq = lapply(model$obsEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsVar = lapply(model$obsVar, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  }
  
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
  
  #Initialize C++ file
  full_modelname = paste(model$modelname,".cpp",sep="")
  fileconn = file(full_modelname)
  
  txt = "#include <TMB.hpp>"
  txt = append(txt, "using namespace density;")
  writeLines(txt,full_modelname)
  
  # Find time-dep inputs
  timedepInput = sort(unique(unlist(lapply(sdeEq,all.vars))))
  timedepInput = timedepInput[!timedepInput %in% diff.processes]
  timedepInput = timedepInput[!timedepInput %in% paste("d",state,sep="")]
  timedepInput = timedepInput[!timedepInput %in% names(c(data$pars, data$constants))]
  
  # Constructing augmented matrices for 1-step solution
  A    = matrix(0,nrow=n,ncol=n)
  negA = matrix(0,nrow=n,ncol=n)
  G    = matrix(0,nrow=n,ncol=n)
  GGT  = matrix(0,nrow=n,ncol=n)
  B    = rep(0,n)
  VarMat.vars = c()
  MeanMat.vars = c()
  for(i in 1:n){
    # Elements in B vector can be found by setting all states equal to zero in the drift term
    mylist = setNames(as.vector(rep(0,n),mode="list"),state)
    cur.expr = Simplify(do.call(substitute,list(diff.terms[[i]]$dt,mylist)))
    MeanMat.vars = c(MeanMat.vars, all.vars(cur.expr))
    B[i] = paste(deparse(hat2pow(cur.expr)),collapse="")
    for(j in 1:n){
      # Elements in A matrix can be found by taking remaining terms after individual state derivatives
      cur.expr = sapply(state, function(x) D(diff.terms[[i]]$dt,x))[[j]]
      VarMat.vars = c(VarMat.vars , all.vars(cur.expr))
      MeanMat.vars = c(MeanMat.vars, all.vars(cur.expr))
      A[i,j] = paste(deparse(hat2pow(cur.expr)),collapse="")
      negA[i,j] = paste(deparse(hat2pow(Simplify(substitute( -a, list(a=cur.expr))))),collapse="")
    }
    # Elements in G matrix are always just diffusion constants in the diagonal
    cur.expr = diff.terms[[i]][[i+1]]
    VarMat.vars = c(VarMat.vars , all.vars(cur.expr))
    G[i,i] = paste(deparse(hat2pow(cur.expr)),collapse="")
    GGT[i,i] = paste(G[i,i],G[i,i],sep="*")
  }
  #Mean
  MeanAug = matrix(0,nrow=n+1,ncol=n+1)
  MeanAug[1:n,1:n] = A
  MeanAug[1:n,n+1] = B
  #Var
  VarAug = matrix(0,nrow=2*n,ncol=2*n)
  VarAug[1:n,1:n] = negA
  VarAug[1:n,(n+1):(2*n)] = GGT
  VarAug[(n+1):(2*n),(n+1):(2*n)] = t(A)
  
  VarMat.vars0 = sort(unique(VarMat.vars))
  VarMat.vars1 = paste("Type", VarMat.vars0, collapse=", ")
  VarMat.vars2 = VarMat.vars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      VarMat.vars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=VarMat.vars2)
    }
  }
  
  MeanMat.vars0 = sort(unique(MeanMat.vars))
  MeanMat.vars1= paste("Type", MeanMat.vars0, collapse=", ")
  MeanMat.vars2 = MeanMat.vars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      MeanMat.vars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=MeanMat.vars2)
    }
  }
  
  temptxt = sprintf("template<class Type>\nmatrix<Type> VarMat(%s){",VarMat.vars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> G(%i,%i);",2*n,2*n))
  for(i in 1:(2*n)){
    for(j in 1:(2*n)){
      temptxt = append( temptxt , sprintf("\tG(%i,%i) = %s;",i-1,j-1,VarAug[i,j]))
    }
  }
  temptxt = append(temptxt, "\treturn G;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  temptxt = sprintf("template<class Type>\nmatrix<Type> MeanMat(%s){",MeanMat.vars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> A(%i,%i);",n+1,n+1))
  for(i in 1:(n+1)){
    for(j in 1:(n+1)){
      temptxt = append( temptxt , sprintf("\tA(%i,%i) = %s;",i-1,j-1,MeanAug[i,j]))
    }
  }
  temptxt = append(temptxt, "\treturn A;\n}")
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
  temptxt = sprintf("\tvector<Type> Xi(%i);",n)
  temptxt = append(temptxt, sprintf("\tvector<Type> Xip1(%i);",n))
  temptxt = append(temptxt, sprintf("\tvector<Type> Z(%i);",n))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> MeanAug0(%i,%i);",n+1,n+1))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> VarAug0(%i,%i);",2*n,2*n))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> MeanAug(%i,%i);",n+1,n+1))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> VarAug(%i,%i);",2*n,2*n))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> Ahat(%i,%i);",n,n))
  temptxt = append(temptxt, sprintf("\tvector<Type> Bhat(%i);",n))
  temptxt = append(temptxt, sprintf("\tmatrix<Type> Qhat(%i,%i);",n,n))
  temptxt = append(temptxt, sprintf("\tType dt;"))
  
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  flag = length(timedepInput[!timedepInput %in% state])<1 & (diff(data$t)[1] - 1e-10) < diff(data$t) && diff(data$t) < (diff(data$t)[1] + 1e-10) 
  if(flag){
    txt = append(txt, sprintf("\tdt = t(1)-t(0);"))
    txt = append(txt, sprintf("\tMeanAug0 = MeanMat(%s)*dt;",paste(MeanMat.vars2,collapse=", ")))
    txt = append(txt, sprintf("\tVarAug0 = VarMat(%s)*dt;",paste(VarMat.vars2,collapse=", ")))
    txt = append(txt, sprintf("\tMeanAug = expm(MeanAug0);"))
    txt = append(txt, sprintf("\tVarAug = expm(VarAug0);"))
    txt = append(txt, sprintf("\tAhat = MeanAug.block(0,0,%i,%i);",n,n))
    txt = append(txt, sprintf("\tBhat = MeanAug.col(%i).head(%i);",n,n))
    txt = append(txt, sprintf("\tQhat = VarAug.block(%i,%i,%i,%i).transpose() * VarAug.block(0,%i,%i,%i);",n,n,n,n,n,n,n,n) )
  }
  
  # dt for loop
  txt = append(txt,"\n\tfor(int i=0;i<t.size()-1;i++){" )
  writeLines(txt,full_modelname)
  if(!flag){
    txt = append(txt, sprintf("\t\tdt = t(i+1) - t(i);"))
    txt = append(txt, sprintf("\t\tMeanAug0 = MeanMat(%s)*dt;",paste(MeanMat.vars2,collapse=", ")))
    txt = append(txt, sprintf("\t\tVarAug0 = VarMat(%s)*dt;",paste(VarMat.vars2,collapse=", ")))
    txt = append(txt, sprintf("\t\tMeanAug = expm(MeanAug0);"))
    txt = append(txt, sprintf("\t\tVarAug = expm(VarAug0);"))
    txt = append(txt, sprintf("\t\tAhat = MeanAug.block(0,0,%i,%i);",n,n))
    txt = append(txt, sprintf("\t\tBhat = MeanAug.col(%i).head(%i);",n,n))
    txt = append(txt, sprintf("\t\tQhat = VarAug.block(%i,%i,%i,%i).transpose() * VarAug.block(0,%i,%i,%i);",n,n,n,n,n,n,n,n) )
  }
  txt = append(txt, sprintf("\t\tXi << %s;",paste(state,"(i)",collapse=", ",sep="")))
  txt = append(txt, sprintf("\t\tXip1 << %s;",paste(state,"(i+1)",collapse=", ",sep="")))
  txt = append(txt, sprintf("\t\tZ = Xip1 - (Ahat * Xi + Bhat);"))
  txt = append(txt, "\t\tnll += MVNORM(Qhat)(Z);\n\t}")
  writeLines(txt,full_modelname)
  
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
