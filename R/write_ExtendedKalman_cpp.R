write_ExtendedKalman_cpp = function(model,data) {
  
  
  # Substitute algebraic expressions
  for(i in 1:length(model$algeqs)){
    curlist = list()
    curlist[[names(model$algeqs)[i]]] = model$algeqs[[i]][[1]]
    model$sdeEq = lapply(model$sdeEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsEq = lapply(model$obsEq, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
    model$obsVar = lapply(model$obsVar, function(x) as.expression(do.call("substitute",list(x[[1]],curlist))))
  }
  
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
  full_modelname = paste(model$modelname,".cpp",sep="")
  fileconn = file(full_modelname)
  
  txt = "#include <TMB.hpp>"
  txt = append(txt, "#include <cmath>")
  txt = append(txt, "using namespace density;")
  txt = append(txt, "template <class Type>
vector<Type> isNA(vector<Type> x){
  int n = x.size();
	vector<Type> NA_boolvec(n);
	NA_boolvec.fill(1.0);
	for(int i=0;i<n;i++){
		if(R_IsNA(asDouble(x(i)))){
			NA_boolvec(i) = 0.0;
		}
	}
	return NA_boolvec;
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
    fvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),fvars2)
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
    dfdxvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),dfdxvars2)
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
    gvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),gvars2)
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
    hvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),hvars2)
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
    dhdxvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),dhdxvars2)
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
  txt = append(txt,"\tint N;")
  txt = append(txt, "\tint s;")
  writeLines(txt,full_modelname)
  
  # Likelihood
  txt = append(txt,"\tType nll = 0;\n")
  writeLines(txt,full_modelname)
  
  txt = append(txt,"\tvector<Type> x0 = X0;\n\tmatrix<Type> p0 = P0;\n\tvector<Type> x1 = x0;\n\tmatrix<Type> p1 = p0;")
  writeLines(txt,full_modelname)
  
  txt = append(txt, "\tint kkk = CppAD::Integer((diff(t)/dt).sum() + 1);")
  writeLines(txt,full_modelname)
  
  k = sum(diff(data$t)/data$dt)+1
  txt = append(txt, 
               sprintf("\tvector<vector<Type>> xPrior(t.size());
\tvector<vector<Type>> xPost(t.size());
\tvector<matrix<Type>> pPrior(t.size());
\tvector<matrix<Type>> pPost(t.size());
\tvector<matrix<Type>> Erep(t.size());
\tvector<matrix<Type>> Crep(t.size());
\tvector<matrix<Type>> C2rep(t.size());
\tvector<matrix<Type>> V2rep(t.size());
\tvector<matrix<Type>> Rrep(t.size());
\tvector<matrix<Type>> Rirep(t.size());
\tvector<matrix<Type>> Krep(t.size());
\tvector<vector<Type>> e0rep(t.size());
\tvector<vector<Type>> erep(t.size());
\tvector<vector<Type>> xPrior_all(kkk);
\tvector<matrix<Type>> pPrior_all(kkk);
\tvector<vector<Type>> F_all(kkk);
\tvector<matrix<Type>> A_all(kkk);
\tvector<matrix<Type>> G_all(kkk);
\txPrior(0) = x0;
\tpPrior(0) = p0;									
\txPost(0) = x0;
\tpPost(0) = p0;
\tmatrix<Type> C,R,K,E,C2,V2,Ri,A,G;	
\tvector<Type> e0,e,y(%s),NAs,F,ynew;
\tmatrix<Type> I(%s,%s);
\tI.setIdentity();
\tmatrix<Type> V(%s,%s);
\tV.setZero();",m,n,n,m,m))
  writeLines(txt,full_modelname)
  
  # Observation variance diagonal
  varobs = paste(unlist(lapply(model$obsVar,function(x) deparse(x[[1]]))),collapse=", ")
  txt = append(txt, sprintf("\tV.diagonal() << %s;",varobs))
  writeLines(txt,full_modelname)
  
  # integration for-loop
  txt = append(txt, "\n\tfor(int i=0;i<t.size()-1;i++){
\t\tN = CppAD::Integer((t(i+1)-t(i))/dt);
\t\tfor(int j=0;j<N;j++){")
  txt = append(txt, sprintf("\t\t\tA  = dfdx(%s);",paste(dfdxvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\tF  = f(%s);",paste(fvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\tG  = g(%s);",paste(gvars2,collapse=", ")) )
  txt = append(txt,"\t\t\tx1 = x0 + F * dt;
			p1 = p0 + ( A*p0 + p0*A.transpose() + G*G.transpose() ) * dt;
			x0 = x1;
			p0 = p1;
			xPrior_all(i*N+j) = x0;
			pPrior_all(i*N+j) = p0;
			F_all(i*N+j) = F;
			A_all(i*N+j) = A;
			G_all(i*N+j) = G;
		}
		xPrior(i+1) = x0;
		pPrior(i+1) = p0;")
  
  obsvars0 = unlist(lapply(obsEq,function(x) deparse(x[[1]][[2]])))
  obsvars2 = paste(obsvars0,"(i)",sep="")
  
  txt = append(txt, sprintf("\t\ty << %s;", paste(obsvars2,collapse=", ")))
  txt = append(txt, sprintf("\t\tNAs = isNA(y);
    s = CppAD::Integer(sum(NAs));
		/*if all are NA we skip, otherwise update (reduce dimensions with E matrix)*/
		if( s > 0 ){
		  ynew = removeNAs(y,NAs);
			matrix<Type> E0(s,%s);
			E 	= constructE(E0,NAs);",m))
  txt = append(txt, sprintf("\t\t\te0  = ynew - h(%s);",paste(hvars2,collapse=", ")))
  txt = append(txt, "\t\t\te 	= E*e0;")
  txt = append(txt, sprintf("\t\t\tC  	= dhdx(%s);",paste(dhdxvars2,collapse=", ")))
  txt = append(txt,"\t\t\tC2 	= E*C;
			V2  = E*V*E.transpose();
			R 	= C2*p0*C2.transpose() + V2;
			Ri  = R.inverse();
			K 	= p0 * C2.transpose() * Ri;
			x0  = x0 + K*e;
			p0  = (I-K*C2)*p0*(I-K*C2).transpose() + K*V2*K.transpose();
			nll += 0.5*atomic::logdet(R) + 0.5*(e*(Ri*e)).sum() + 0.5*log(2*M_PI)*asDouble(s);
			Erep(i) = E;
			e0rep(i) = e0;
			erep(i) = e;
			Crep(i) = C;
			C2rep(i) = C2;
			V2rep(i) = V2;
			Rrep(i) = R;
			Rirep(i) = Ri;
			Krep(i) = K;
		}
		xPost(i+1) = x0;
		pPost(i+1) = p0;
  }
REPORT(xPrior);
REPORT(xPost);
REPORT(pPrior);
REPORT(pPost);
REPORT(xPrior_all);
REPORT(pPrior_all);
REPORT(F_all);
REPORT(A_all);
REPORT(G_all);
REPORT(Erep);
REPORT(Crep);
REPORT(C2rep);
REPORT(Rrep);
REPORT(Rirep);
REPORT(V2rep);
REPORT(Krep);
REPORT(e0rep);
REPORT(erep);
return nll;
}")
  writeLines(txt,full_modelname)
  
  # Close file connection
  close(fileconn)
  
}
