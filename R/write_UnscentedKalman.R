write_UnscentedKalman = function(model,data) {
  
  #Get info
  model_info = getModelInfo(model,data)
  names_model_info = names(model_info)
  for(i in 1:length(model_info)){
    assign(names_model_info[i],model_info[[i]])
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
  
  ##################################################
  # Various helper functions
  ##################################################
  
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
  
  txt = append(txt,"template<class Type>
vector<Type> removeNAs(int s, vector<Type> y,vector<Type> NAs){
  int n = y.size();
  int count = 0;
	vector<Type> ans(s);
	for(int i=0;i<n;i++){
		if(NAs(i) > 0.5){
			ans(count) = y(i);
			count += 1;
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
  
  temptxt = sprintf("template<class Type>
matrix<Type> construct_Xsp(vector<Type> x0, matrix<Type> s0){
  int n = %s;
  int nn = %s;
  matrix<Type> Xsp(n,nn);
  vector<Type> Si;
  Xsp.col(0) = x0;
  for(int i=1; i<n+1; i++){
    Si = s0.col(i-1);
    Xsp.col(i) = x0 + sqrt(1+n) * Si;
    Xsp.col(i+n) = x0 - sqrt(1+n) * Si;
  }
  return Xsp;
}",n,2*n+1)
  txt = append(txt,temptxt)
  
  txt = append(txt, "template<class Type>
matrix<Type> Phi(matrix<Type> M){
  matrix<Type> K(M.col(0).size(),M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/Type(2.0);
  return K;
}")
  
  writeLines(txt,full_modelname)
  
  # Construct drift f, diffusion g, observation h
  
  ##################################################
  # Construct drift function
  ##################################################
  
  fvars0 = sort(unique(unlist(lapply(1:n,function(i) all.vars(diff.terms[[i]]$dt)))))
  fvars0_sigma = fvars0[!(fvars0 %in% state)]
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars1_sigma = paste("Type", fvars0_sigma, collapse=", ")
  fvars2 = fvars0
  fvars2_sigma = fvars0
  fvars2_sigma2 = fvars0_sigma
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      fvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=fvars2)
      fvars2_sigma2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=fvars2_sigma2)
    }
  }
  for(i in 1:n){
    fvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),fvars2)
    fvars2_sigma = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),fvars2_sigma)
  }
  if(length(fvars2)<1){
    fvars1 = "Type a"
    fvars2 = "Type(0.0)"
    fvars1_sigma = "Type a"
    fvars2_sigma = "Type(0.0)"
    fvars2_sigma2 = "Type(0.0)"
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
  # Construct sigma points drift function
  ##################################################
  
  temptxt = sprintf("template<class Type>
matrix<Type> construct_F(matrix<Type> Xsp, %s){
  int n = %s;
  int nn = %s;
  matrix<Type> F(n,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    F.col(i) = f(%s);
  }
  return F;
}",fvars1_sigma,n,2*n+1,paste(fvars2_sigma,collapse=", "))
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Construct diffusion function
  ##################################################
  
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
    # gvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),gvars2)
    gvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__Xsp0(0,%s)",i-1),gvars2)
  }
  
  g = matrix(0,nrow=n,ncol=ng)
  for(i in 1:n){
    for(j in 1:ng){
      term = paste(deparse(hat2pow(diff.terms[[i]][[j+1]])),collapse = "")
      g[i,j] = term
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> g(%s){",gvars1)
  temptxt = append(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",n,ng))
  for(i in 1:n){
    for(j in 1:ng){
      temptxt = append(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,g[i,j]))
    }
  }
  temptxt = append(temptxt, "\treturn ans;\n}")
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Construct observation function
  ##################################################
  
  hvars0 = sort(unique(unlist(lapply(lapply(obsEq,function(x) x[[1]][[3]]), all.vars))))
  hvars0_sigma = hvars0[!(hvars0 %in% state)]
  hvars1 = paste("Type", hvars0, collapse=", ")
  hvars1_sigma = paste("Type", hvars0_sigma, collapse=", ")
  hvars2 = hvars0
  hvars2_sigma = hvars0
  hvars2_sigma2 = hvars0_sigma
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      hvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=hvars2)
      hvars2_sigma2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=hvars2_sigma2)
    }
  }
  for(i in 1:n){
    hvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("__x0(%s)",i-1),hvars2)
    hvars2_sigma = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0(%s)",i-1),hvars2_sigma)
  }
  if(length(hvars0_sigma)<1){
    hvars1_sigma = "Type a"
    hvars2_sigma2 = "Type(0.0)"
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
  # Construct observation sigma points function
  ##################################################
  
  temptxt = sprintf("template<class Type>
matrix<Type> construct_H(matrix<Type> Xsp, %s){
  int nn = %s;
  int m = %s;
  matrix<Type> H(m,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    H.col(i) = h(%s);
  }
  return H;
}",hvars1_sigma,2*n+1,m, paste(hvars2_sigma,collapse=", "))
  txt = append(txt, temptxt)
  writeLines(txt,full_modelname)
  
  ##################################################
  # Write TMB objective function
  ##################################################
  
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
  
  txt = append(txt,"\t/* DEFINE VARIABLES*/")
  
  txt = append(txt, sprintf("\tint __n = %s;",n))
  txt = append(txt, sprintf("\tint __nn = %s;",2*n+1))
  txt = append(txt, sprintf("\tint __m = %s;",m))
  txt = append(txt,"\tint __N,__s;")
  txt = append(txt,"\tType __nll = 0;\n")
  txt = append(txt, "\tvector<Type> __y(__m);")
  txt = append(txt,"\tvector<Type> __NAs,__ynew,__e,__X1;")
  txt = append(txt, "\tmatrix<Type> __Xsp0,__S0,__G, __F, __S0Inv, __M, __S0PhiM, __Frhs1(__n,__nn), __Frhs0, __Frhs;")
  txt = append(txt, "\tmatrix<Type> __E, __H, __Syy, __SyyInv, __Sxy, __K, __P1;")
  txt = append(txt, "\tvector<vector<Type>> xPost(t.size());")
  txt = append(txt, "\tvector<matrix<Type>> pPost(t.size());")
  txt = append(txt, "\tvector<vector<Type>> OneStepErrors(t.size());")
  txt = append(txt, "\tvector<matrix<Type>> OneStepErrorsCovariance(t.size());")
  txt = append(txt, "\txPost(0) = X0;")
  txt = append(txt, "\tpPost(0) = P0;\n")
  writeLines(txt,full_modelname)
  
  # Observation variance matrix
  txt = append(txt,"\tmatrix<Type> __V(__m,__m);")
  txt = append(txt,"\t__V.setZero();")
  varobs = paste(unlist(lapply(model$obsVar,function(x) deparse(x[[1]]))),collapse=", ")
  txt = append(txt, sprintf("\t__V.diagonal() << %s;",varobs))
  writeLines(txt,full_modelname)
  
  txt = append(txt ,"\n\t/*CONSTRUCT WEIGHTS*/
  Type __lambda = 3-__n, __weights;
  vector<Type> __Wm(__nn);
  matrix<Type> __Wmm(__nn,__nn), __Wm_diag(__nn,__nn), __I(__nn,__nn), __W;
  __weights = Type(1.0)/(Type(2.0)*(__lambda+__n));
  __Wm.fill(__weights);
  __Wm(0) = __lambda/(__lambda+__n);
  for(int i=0; i<__nn; i++){
    __Wmm.col(i) = __Wm;
  }
  __Wm_diag.setZero(); 
  __Wm_diag.diagonal() = __Wm;
  __I.setIdentity();
  __W = (__I-__Wmm) * __Wm_diag * (__I-__Wmm).transpose();\n")
  writeLines(txt,full_modelname)
  
  # integration for-loop
  txt = append(txt,"\n\t/* MAIN FOR-LOOP - FILTERING/UPDATE*/")
  txt = append(txt, "\t __S0 = P0.llt().matrixL();")
  txt = append(txt, "\t __Xsp0 = construct_Xsp(X0,__S0);")
  txt = append(txt, "\t for(int i=0;i<t.size()-1;i++){
\t\t __N = CppAD::Integer((t(i+1)-t(i))/dt);
\t\t for(int j=0;j<__N;j++){")
  txt = append(txt, sprintf("\t\t\t __G  = g(%s);",paste(gvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\t __F  = construct_F(__Xsp0,%s);",paste(fvars2_sigma2,collapse=", ")))
  txt = append(txt, "\t\t\t __S0Inv = __S0.inverse();")
  txt = append(txt, "\t\t\t __M = __S0Inv * (__Xsp0 * __W * __F.transpose() + __F * __W * __Xsp0.transpose() + __G*__G.transpose()) * __S0Inv.transpose();")
  txt = append(txt, "\t\t\t __S0PhiM = __S0 * Phi(__M);")
  txt = append(txt, "\t\t\t __Frhs1.block(0,1,__n,__n) = __S0PhiM;")
  txt = append(txt, "\t\t\t __Frhs1.block(0,__n+1,__n,__n) = -__S0PhiM;")
  txt = append(txt, "\t\t\t __Frhs0 = (__F*__Wm).replicate(1,__nn);")
  txt = append(txt, "\t\t\t __Frhs = __Frhs0 + sqrt(3.0) * __Frhs1;")
  txt = append(txt, "\t\t\t __Xsp0 += __Frhs * dt;")
  txt = append(txt, "\t\t\t __S0 = ((__Xsp0 - __Xsp0.col(0).replicate(1,__nn))/sqrt(Type(3.0))).block(0,1,__n,__n);")
  txt = append(txt, "\t\t\t};")
  writeLines(txt,full_modelname)
  
  obsvars0 = unlist(lapply(obsEq,function(x) deparse(x[[1]][[2]])))
  obsvars2 = paste(obsvars0,"(i+1)",sep="")
  
  txt = append(txt, "\t\t__P1 = __S0 * __S0.transpose();")
  txt = append(txt, sprintf("\t\t__y << %s;", paste(obsvars2,collapse=", ")))
  txt = append(txt, "\t\t__NAs = isNA(__y);")
  txt = append(txt, "\t\t__s = CppAD::Integer(sum(__NAs));")
  txt = append(txt, "\t\tif( __s > 0 ){")
  
  txt = append(txt, "\t\t\t __ynew = removeNAs(__s,__y,__NAs);")
  txt = append(txt, "\t\t\t matrix<Type> __E0(__s,__m);")
  txt = append(txt, "\t\t\t __E 	= constructE(__E0,__NAs);")
  txt = append(txt, sprintf("\t\t\t __H = construct_H(__Xsp0,%s);",hvars2_sigma2))
  txt = append(txt, "\t\t\t __e  = __ynew - __E * (__H * __Wm);")
  txt = append(txt, "\t\t\t __Syy  = __E * (__H * __W * __H.transpose() + __V) * __E.transpose();")
  txt = append(txt, "\t\t\t __Sxy  = __Xsp0 * __W * __H.transpose() * __E.transpose();")
  txt = append(txt, "\t\t\t __SyyInv  = __Syy.inverse();")
  txt = append(txt, "\t\t\t __K = __Sxy * __SyyInv;")
  txt = append(txt, "\t\t\t __X1 = __Xsp0 * __Wm + __K * __e;")
  txt = append(txt, "\t\t\t __P1 = __S0 * __S0.transpose() - __K * __Syy * __K.transpose();")
  txt = append(txt, "\t\t\t __S0 = __P1.llt().matrixL();")
  txt = append(txt, "\t\t\t __Xsp0 = construct_Xsp(__X1,__S0);")
  txt = append(txt, "\t\t\t __nll += 0.5*atomic::logdet(__Syy) + 0.5*(__e*(__SyyInv*__e)).sum() + 0.5*log(2*M_PI)*asDouble(__s);")
  txt = append(txt, "\t\t\t OneStepErrors(i) = __e;")
  txt = append(txt, "\t\t\t OneStepErrorsCovariance(i) = __Syy;")
  txt = append(txt, "\t\t};")
  txt = append(txt, "\t\t xPost(i+1) = __Xsp0.col(0);")
  txt = append(txt, "\t\t pPost(i+1) = __P1;")
  txt = append(txt, "\t};")
  txt = append(txt,"\tREPORT(xPost);")
  txt = append(txt,"\tREPORT(pPost);")
  txt = append(txt,"\tREPORT(OneStepErrors);")
  txt = append(txt,"\tREPORT(OneStepErrorsCovariance);")
  txt = append(txt, "\treturn __nll;")
  txt = append(txt, "};")
  writeLines(txt,full_modelname)
  
  # Close file connection
  close(fileconn)
}
