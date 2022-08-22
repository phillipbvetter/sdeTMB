write_ExtendedKalman_cpp = function(model,data) {
  
  
  # Substitute algebraic expressions
  # for(i in 1:length(model$algeqs)){
  #   curlist = list()
  #   curlist[[names(model$algeqs)[i]]] = model$algeqs[[i]][[1]]
  #   model$sdeEq = lapply(model$sdeEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  #   model$obsEq = lapply(model$obsEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  #   model$obsVar = lapply(model$obsVar, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  # }
  
  # Extract
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
  timedepInput = timedepInput[!timedepInput %in% state]
  
  #Initialize C++ file
  full_modelname = paste(model$modelname2,".cpp",sep="")
  fileconn = file(full_modelname)
  
  txt = "#include <TMB.hpp>"
  txt = append(txt, "#include <cmath>")
  txt = append(txt, "using namespace density;")
  txt = append(txt, "template <class Type>
vector<Type> isNA(vector<Type> x){
  int n = x.size();
	vector<Type> ans(n);
	ans.fill(1.0);
	for(int i=0;i<n;i++){
		if(R_IsNA(asDouble(x(i)))){
			ans(i) = 0.0;
		}
	}
	return ans;
}")
  txt = append(txt, "template<class Type>
vector<Type> removeNAs(vector<Type> y,vector<Type> NAs){
  int n = y.size();
	vector<Type> ans = y;
	for(int i=0;i<n;i++){
		if(NAs(i) < 0.5){
			ans(i) = 0.0;
		} else {
			ans(i) = y(i);
		}
	}
  return ans;
}")
  txt = append(txt, "template <class Type>
matrix<Type> constructE(matrix<Type> E,vector<Type> p){
	matrix<Type> ans = E;
	ans.setZero();
	int s = E.row(0).size();
	int j=0;
	for(int i=0;i<s;i++){
		if(p(i) == 1){ /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
			ans(j,i) = 1.0;
			j += 1;
		}
	}
	return ans;
}")
  writeLines(txt,full_modelname)
  
  # Construct drift f, jacobian df/dx, diffusion g, observation h and jacobian dh/dx
  ##################################################
  # Construct drift function
  
  fvars0 = sort(unique(unlist(lapply(lapply(rhs,function(x) D(x,"dt")),all.vars))))
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars2 = fvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      fvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=fvars2)
    }
  }
  for(i in 1:n){
    fvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),fvars2)
  }
  if(length(fvars2)<1){
    fvars1 = "Type a"
    fvars2 = "Type(0.0)"
  }
  
  temptxt = sprintf("template<class Type>\nvector<Type> f(%s){",fvars1)
  temptxt = append(temptxt, sprintf("\tvector<Type> ans(%i);",n))
  for(i in 1:n){
    term = paste(deparse(hat2pow(diff.terms[[i]]$dt)),collapse="")
    temptxt = append(temptxt, sprintf("\tans(%i) = %s;",i-1,term))
  }
  temptxt = append(temptxt, "\treturn ans;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Construct jacobian of drift function
  
  dfdxvars0 = c()
  for(i in 1:n){
    dfdxvars0 = c(dfdxvars0, unlist(lapply(lapply(state,function(x) D(diff.terms[[i]]$dt,x)),all.vars)))
  }
  dfdxvars0 = sort(unique(dfdxvars0))
  dfdxvars1 = paste("Type", dfdxvars0, collapse=", ")
  dfdxvars2 = dfdxvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      dfdxvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=dfdxvars2)
    }
  }
  for(i in 1:n){
    dfdxvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),dfdxvars2)
  }
  if(length(dfdxvars2)<1){
    dfdxvars1 = "Type a"
    dfdxvars2 = "Type(0.0)"
  }
  
  dfdx = matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    terms = lapply(state,function(x) hat2pow(D(diff.terms[[i]]$dt,x)))
    for(j in 1:n){
      dfdx[i,j] = paste(deparse(terms[[j]]),collapse="")
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> dfdx(%s){",dfdxvars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",n,n))
  for(i in 1:n){
    for(j in 1:n){
      temptxt = append(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,dfdx[i,j]))
    }
  }
  temptxt = append(temptxt, "\treturn ans;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Construct diffusion function
  gvars0 = c()
  for(i in 1:n){
    gvars0 = c(gvars0, unlist(lapply(diff.terms[[i]][-1],all.vars)))
  }
  gvars0 = sort(unique(gvars0))
  gvars1 = paste("Type", gvars0, collapse=", ")
  gvars2 = gvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      gvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=gvars2)
    }
  }
  for(i in 1:n){
    gvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),gvars2)
  }
  
  g = matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    for(j in 1:n){
      term = paste(deparse(hat2pow(diff.terms[[i]][[j+1]])),collapse = "")
      g[i,j] = term
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> g(%s){",gvars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",n,n))
  for(i in 1:n){
    for(j in 1:n){
      temptxt = append(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,g[i,j]))
    }
  }
  temptxt = append(temptxt, "\treturn ans;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Construct observation function
  
  hvars0 = unlist(lapply(lapply(obsEq,function(x) x[[1]][[3]]), all.vars))
  hvars0 = sort(unique(hvars0))
  hvars1 = paste("Type", hvars0, collapse=", ")
  hvars2 = hvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      hvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=hvars2)
    }
  }
  for(i in 1:n){
    hvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),hvars2)
  }
  
  h = c()
  for(i in 1:m){
    term = paste(deparse(hat2pow(obsEq[[i]][[1]][[3]])),collapse = "")
    h[i] = term
  }
  temptxt = sprintf("template<class Type>\nvector<Type> h(%s){",hvars1)
  temptxt = append(temptxt, sprintf("\tvector<Type> ans(%i);",m))
  for(i in 1:m){
    temptxt = append(temptxt, sprintf("\tans(%s) = %s;",i-1,h[i]))
  }
  temptxt = append(temptxt, "\treturn ans;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Construct jacobian of observation function
  
  dhdxvars0 = c()
  for(i in 1:m){
    dhdxvars0 = c(dhdxvars0 , unlist(lapply(lapply(state, function(x) D(obsEq[[i]][[1]][[3]],x)),all.vars)))
  }
  dhdxvars0 = sort(unique(dhdxvars0))
  dhdxvars1 = paste("Type", dhdxvars0, collapse=", ")
  dhdxvars2 = dhdxvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      dhdxvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=dhdxvars2)
    }
  }
  for(i in 1:n){
    dhdxvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),dhdxvars2)
  }
  if(length(dhdxvars2)<1){
    dhdxvars1 = "Type a"
    dhdxvars2 = "Type(0.0)"
  }
  
  dhdx = matrix(NA,nrow=m,ncol=n)
  for(i in 1:m){
    terms = lapply(state, function(x) hat2pow(D(obsEq[[i]][[1]][[3]],x)))
    for(j in 1:n){
      dhdx[i,j] = paste(deparse(terms[[j]]),collapse = "")
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> dhdx(%s){",dhdxvars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",m,n))
  for(i in 1:m){
    for(j in 1:n){
      temptxt = append(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,dhdx[i,j]))
    }
  }
  temptxt = append(temptxt, "\treturn ans;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  
  txt = append(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")
  
  NOTdataNames = c("constants","pars",state,"X0","P0","dt")
  dataNames = names(data)[!(names(data) %in% NOTdataNames)]
  # Data and Parameters
  for(i in 1:length(dataNames)){
    nam = dataNames[i]
    txt = append(txt, sprintf("\tDATA_VECTOR(%s);",nam),length(txt))
  }
  writeLines(txt,full_modelname)
  
  # Parameters
  for(i in 1:length(data$pars)){
    nam = names(data$pars)[i]
    txt = append(txt, sprintf("\tPARAMETER(%s);",nam))
  }
  writeLines(txt,full_modelname)
  
  # Various matrices and vectors
  txt = append(txt, "\tDATA_VECTOR(X0);")
  txt = append(txt, "\tDATA_MATRIX(P0);")
  txt = append(txt, "\tDATA_SCALAR(dt);")
  writeLines(txt,full_modelname)  

  # Constants
  for(i in 1:length(data$constants)){
    nam = names(data$constants)[i]
    if(length(data$constants[[i]])>1){
      txt = append( txt , sprintf("\tDATA_VECTOR(%s);",nam))
    } else {
      txt = append( txt , sprintf("\tDATA_SCALAR(%s);",nam))
    }
  }
  txt = append(txt,"\tint __N;")
  txt = append(txt, "\tint __s;")
  writeLines(txt,full_modelname)
  
  # Likelihood
  txt = append(txt,"\tType __nll = 0;\n")
  writeLines(txt,full_modelname)
  
  txt = append(txt,"\tvector<Type> __x0 = X0;\n\tmatrix<Type> __p0 = P0;\n\tvector<Type> __x1 = __x0;\n\tmatrix<Type> __p1 = __p0;")
  writeLines(txt,full_modelname)
  
  txt = append(txt, "\tint __k = CppAD::Integer((diff(t)/dt).sum() + 1);")
  writeLines(txt,full_modelname)
  
  k = sum(diff(data$t)/data$dt)+1
  txt = append(txt, 
               sprintf("
\tvector<vector<Type>> __xPrior(t.size());
\tvector<vector<Type>> __xPost(t.size());
\tvector<vector<Type>> __xPriorPost(2*t.size());
\tvector<matrix<Type>> __pPrior(t.size());
\tvector<matrix<Type>> __pPost(t.size());
\tvector<matrix<Type>> __pPriorPost(2*t.size());
\tvector<vector<Type>> __xPrior_all(__k);
\tvector<matrix<Type>> __pPrior_all(__k);
\tvector<vector<Type>> __F_all(__k);
\tvector<matrix<Type>> __A_all(__k);
\tvector<matrix<Type>> __G_all(__k);
\t__xPrior(0) = __x0;
\t__pPrior(0) = __p0;									
\t__xPost(0) = __x0;
\t__pPost(0) = __p0;
\t__xPriorPost(0) = __x0;
\t__xPriorPost(1) = __x0;
\t__pPriorPost(0) = __p0;
\t__pPriorPost(1) = __p0;

\tmatrix<Type> __C,__R,__K,__E,__C2,__V2,__Ri,__A,__G;	
\tvector<Type> __e0,__e,__y(%s),__NAs,__F,__ynew;
\tmatrix<Type> __I(%s,%s);
\t__I.setIdentity();
\tmatrix<Type> __V(%s,%s);
\t__V.setZero();",m,n,n,m,m))
  writeLines(txt,full_modelname)
  
  # Observation variance diagonal
  varobs = paste(unlist(lapply(model$obsVar,function(x) deparse(x[[1]]))),collapse=", ")
  txt = append(txt, sprintf("\t__V.diagonal() << %s;",varobs))
  writeLines(txt,full_modelname)
  
  # integration for-loop
  txt = append(txt, "\n\tfor(int i=0;i<t.size()-1;i++){
\t\t__N = CppAD::Integer((t(i+1)-t(i))/dt);
\t\tfor(int j=0;j<__N;j++){")
  txt = append(txt, sprintf("\t\t\t__A  = dfdx(%s);",paste(dfdxvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\t__F  = f(%s);",paste(fvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\t__G  = g(%s);",paste(gvars2,collapse=", ")) )
  txt = append(txt,"\t\t\t__x1 = __x0 + __F * dt;
			__p1 = __p0 + ( __A*__p0 + __p0*__A.transpose() + __G*__G.transpose() ) * dt;
			__x0 = __x1;
			__p0 = __p1;
			__xPrior_all(i*__N+j) = __x0;
			__pPrior_all(i*__N+j) = __p0;
			__F_all(i*__N+j) = __F;
			__A_all(i*__N+j) = __A;
			__G_all(i*__N+j) = __G;
		}
		__xPrior(i+1) = __x0;
		__pPrior(i+1) = __p0;
		__xPriorPost(2*i+2) = __x0;
		__pPriorPost(2*i+2) = __p0;
		")
  
  obsvars0 = unlist(lapply(obsEq,function(x) deparse(x[[1]][[2]])))
  obsvars2 = paste(obsvars0,"(i)",sep="")
  
  txt = append(txt, sprintf("\t\t__y << %s;", paste(obsvars2,collapse=", ")))
  txt = append(txt, sprintf("\t\t__NAs = isNA(__y);
    __s = CppAD::Integer(sum(__NAs));
		/*if all are NA we skip, otherwise update (reduce dimensions with E matrix)*/
		if( __s > 0 ){
		  __ynew = removeNAs(__y,__NAs);
			matrix<Type> __E0(__s,%s);
			__E 	= constructE(__E0,__NAs);",m))
  txt = append(txt, sprintf("\t\t\t__e0  = __ynew - h(%s);",paste(hvars2,collapse=", ")))
  txt = append(txt, "\t\t\t__e 	= __E*__e0;")
  txt = append(txt, sprintf("\t\t\t__C  	= dhdx(%s);",paste(dhdxvars2,collapse=", ")))
  txt = append(txt,"\t\t\t__C2 	= __E*__C;
			__V2  = __E*__V*__E.transpose();
			__R 	= __C2*__p0*__C2.transpose() + __V2;
			__Ri  = __R.inverse();
			__K 	= __p0 * __C2.transpose() * __Ri;
			__x0  = __x0 + __K*__e;
			__p0  = (__I-__K*__C2)*__p0*(__I-__K*__C2).transpose() + __K*__V2*__K.transpose();
			__nll += 0.5*atomic::logdet(__R) + 0.5*(__e*(__Ri*__e)).sum() + 0.5*log(2*M_PI)*asDouble(__s);
		}
		__xPost(i+1) = __x0;
		__pPost(i+1) = __p0;
		__xPriorPost(2*i+3) = __x0;
    __pPriorPost(2*i+3) = __p0;
  }")
	txt = append(txt ,"
REPORT(__xPrior);
REPORT(__xPost);
REPORT(__xPriorPost);
REPORT(__xPrior_all);
REPORT(__pPrior);
REPORT(__pPost);
REPORT(__pPriorPost);
REPORT(__pPrior_all);
REPORT(__F_all);
REPORT(__A_all);
REPORT(__G_all);
return __nll;
}")
  writeLines(txt,full_modelname)
  
  # Close file connection
  close(fileconn)
  
}
