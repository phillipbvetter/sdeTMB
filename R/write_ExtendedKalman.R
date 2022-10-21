write_ExtendedKalman = function(model,data,control) {
  
  #Get info
  model_info = getModelInfo(model,data)
  names_model_info = names(model_info)
  for(i in 1:length(model_info)){
    assign(names_model_info[i],model_info[[i]])
  }
  
  # Find time-dep inputs
  timedepInput = sort(unique(unlist(c(lapply(sdeEq,all.vars),lapply(obsEq,all.vars),lapply(obsVar,all.vars)))))
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
  
  txt = append(txt,"template<class Type>
Type lossfunction(Type x, vector<Type> tukeypars, Type huber_c, int lossFunc){
  Type loss;
  if(lossFunc==1){
    Type a = tukeypars(0);
    Type b = tukeypars(1);
    Type c = tukeypars(2);
    Type d = tukeypars(3);
    loss = d * ( (Type(1.0)/(Type(1.0)+exp(-a*(x-b)))) + c );
  } else if (lossFunc==2){
    Type c_squared = pow(huber_c,2);
    loss = c_squared * (sqrt(1 + (x / c_squared)) - 1);
  } else {
    loss = x;
  }
  return(loss);
}")
  writeLines(txt,full_modelname)
  
  # Construct drift f, jacobian df/dx, diffusion g, observation h and jacobian dh/dx
  ##################################################
  # Construct drift function
  
  fvars0 = sort(unique(unlist(lapply(1:n,function(i) all.vars(diff.terms[[i]]$dt)))))
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars2 = fvars0
  if(length(timedepInput)>0){
    for(i in 1:length(timedepInput)){
      fvars2 = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i)",timedepInput[i]), x=fvars2)
    }
  }
  for(i in 1:n){
    fvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0__(%s)",i-1),fvars2)
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
    dfdxvars0 = c(dfdxvars0, unlist(lapply(lapply(state,function(x) Deriv(diff.terms[[i]]$dt,x=x,cache.exp=FALSE)),all.vars)))
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
    dfdxvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0__(%s)",i-1),dfdxvars2)
  }
  if(length(dfdxvars2)<1){
    dfdxvars1 = "Type a"
    dfdxvars2 = "Type(0.0)"
  }
  
  dfdx = matrix(0,nrow=n,ncol=n)
  for(i in 1:n){
    terms = lapply(state,function(x) hat2pow(Deriv(diff.terms[[i]]$dt,x=x,cache.exp=FALSE)))
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
    gvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0__(%s)",i-1),gvars2)
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
    hvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0__(%s)",i-1),hvars2)
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
    dhdxvars0 = c(dhdxvars0 , unlist(lapply(lapply(state, function(x) Deriv(obsEq[[i]][[1]][[3]],x=x,cache.exp=FALSE)),all.vars)))
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
    dhdxvars2 = sub(pattern=sprintf("^%s$",state[i]),replacement=sprintf("x0__(%s)",i-1),dhdxvars2)
  }
  if(length(dhdxvars2)<1){
    dhdxvars1 = "Type a"
    dhdxvars2 = "Type(0.0)"
  }
  
  dhdx = matrix(NA,nrow=m,ncol=n)
  for(i in 1:m){
    terms = lapply(state, function(x) hat2pow(Deriv(obsEq[[i]][[1]][[3]],x=x,cache.exp=FALSE)))
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
  
  NOTdataNames = c("constants","pars",state,"X0","P0","dt","MAPmean","MAPcov","MAPbool",
                   "tukeypars","huber_c","loss_int")
  dataNames = names(data)[!(names(data) %in% NOTdataNames)]
  # Data and Parameters
  for(i in 1:length(dataNames)){
    nam = dataNames[i]
    txt = append(txt, sprintf("\t DATA_VECTOR(%s);",nam),length(txt))
  }
  writeLines(txt,full_modelname)
  
  # Parameters
  for(i in 1:length(data$pars)){
    nam = names(data$pars)[i]
    txt = append(txt, sprintf("\t PARAMETER(%s);",nam))
  }
  writeLines(txt,full_modelname)
  
  # Various matrices and vectors
  txt = append(txt, "\t DATA_VECTOR(X0);")
  txt = append(txt, "\t DATA_MATRIX(P0);")
  txt = append(txt, "\t DATA_SCALAR(dt);")
  txt = append(txt, "\t DATA_VECTOR(tukeypars);")
  txt = append(txt, "\t DATA_SCALAR(huber_c);")
  txt = append(txt, "\t DATA_INTEGER(loss_int);")
  txt = append(txt, "\t DATA_INTEGER(MAPbool);")
  writeLines(txt,full_modelname)  
  
  # Constants
  for(i in 1:length(data$constants)){
    nam = names(data$constants)[i]
    if(length(data$constants[[i]])>1){
      txt = append( txt , sprintf("\t DATA_VECTOR(%s);",nam))
    } else {
      txt = append( txt , sprintf("\t DATA_SCALAR(%s);",nam))
    }
  }
  txt = append(txt, sprintf("\n\t int n__ = %s;",n))
  txt = append(txt, sprintf("\t int m__ = %s;",m))
  txt = append(txt, "\t int k__ = CppAD::Integer((diff(t)/dt).sum() + 1);")
  txt = append(txt,"\t int N__;")
  txt = append(txt, "\t int s__;")
  writeLines(txt,full_modelname)
  
  # Likelihood
  txt = append(txt,"\t Type nll__ = 0;")
  writeLines(txt,full_modelname)
  
  txt = append(txt, "\n\t vector<vector<Type>> xPrior(t.size());")
  txt = append(txt, "\t vector<vector<Type>> xPost(t.size());")
  txt = append(txt, "\t vector<vector<Type>> xPriorPost(2*t.size());")
  txt = append(txt, "\t vector<matrix<Type>> pPrior(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> pPost(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> pPriorPost(2*t.size());")
  txt = append(txt, "\t vector<vector<Type>> xPrior_all(k__);")
  txt = append(txt, "\t vector<matrix<Type>> pPrior_all(k__);")
  txt = append(txt, "\t vector<vector<Type>> F_all(k__);")
  txt = append(txt, "\t vector<matrix<Type>> A_all(k__);")
  txt = append(txt, "\t vector<matrix<Type>> G_all(k__);")
  txt = append(txt, "\t vector<vector<Type>> Innovation(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> InnovationCovariance(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> E_mat(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> C_mat(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> V_mat(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> Ri_mat(t.size());")
  txt = append(txt, "\t vector<matrix<Type>> K_mat(t.size());")
  
  txt = append(txt, "\n\t vector<Type> x0__ = X0;")
  txt = append(txt, "\t matrix<Type> p0__ = P0;")
  txt = append(txt, "\t vector<Type> x1__ = x0__;")
  txt = append(txt, "\t matrix<Type> p1__ = p0__;")
  txt = append(txt, "\t xPrior(0) = x0__;")
  txt = append(txt, "\t pPrior(0) = p0__;")
  txt = append(txt, "\t xPost(0) = x0__;")
  txt = append(txt, "\t pPost(0) = p0__;")
  txt = append(txt, "\t xPriorPost(0) = x0__;")
  txt = append(txt, "\t xPriorPost(1) = x0__;")
  txt = append(txt, "\t pPriorPost(0) = p0__;")
  txt = append(txt, "\t pPriorPost(1) = p0__;")
  
  txt = append(txt, "\n\t matrix<Type> C__,R__,K__,E__,C2__,V2__,Ri__,A__,G__;")
  txt = append(txt, "\t vector<Type> e0__,e__,y__(m__),NAs__,F__,ynew__,h__;")
  txt = append(txt, "\t matrix<Type> I__(n__,n__);")
  txt = append(txt, "\t I__.setIdentity();")
  writeLines(txt,full_modelname)
  
  # Observation variance diagonal
  obsvars = lapply(model$obsVar,function(x) all.vars(x[[1]])) #for index substitution
  obsvars0 = obsvars #keeper
  obsvars1 = lapply(model$obsVar,function(x) deparse(x[[1]])) #actual strings to be substituted in
  if(length(timedepInput)>0){
    for(j in 1:m){
      for(i in 1:length(timedepInput)){
        obsvars[[j]] = sub(pattern=sprintf("^%s$",timedepInput[i]), replacement=sprintf("%s(i+1)",timedepInput[i]), x=obsvars[[j]])
      }
      for(i in 1:length(obsvars[[j]])){
        obsvars1[[j]] = sub(pattern=obsvars0[[j]][i], replacement=obsvars[[j]][i], x=obsvars1[[j]])
      }
    }
  }
  txt = append(txt, "\n\t matrix<Type> V__(m__,m__);")
  txt = append(txt, "\t V__.setZero();")
  writeLines(txt,full_modelname)
  
  # Observation equation
  obs.eqs = paste(unlist(lapply(obsEq,function(x) deparse(x[[1]][[2]]))),"(i+1)",sep="")
  
  # integration for-loop
  txt = append(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = append(txt, "\t //////////// MAIN-FOR LOOP OVER OBSERVATIONS ///////////")
  txt = append(txt, "\t ////////////////////////////////////////////////////////")
  
  txt = append(txt, "\t for(int i=0 ; i<t.size()-1 ; i++){")
  txt = append(txt, "\t\t N__ = CppAD::Integer((t(i+1)-t(i))/dt);")
  
  txt = append(txt, "\n\t\t //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  
  txt = append(txt, "\t\t for(int j=0 ; j<N__ ; j++){")
  txt = append(txt, sprintf("\t\t\t A__  = dfdx(%s);",paste(dfdxvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\t F__  = f(%s);",paste(fvars2,collapse=", ")) )
  txt = append(txt, sprintf("\t\t\t G__  = g(%s);",paste(gvars2,collapse=", ")) )
  txt = append(txt,"\t\t\t x1__ = x0__ + F__ * dt;")
  txt = append(txt, "\t\t\t p1__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt;")
  txt = append(txt, "\t\t\t x0__ = x1__;")
  txt = append(txt, "\t\t\t p0__ = p1__;")
  txt = append(txt, "\n\t\t\t xPrior_all(i*N__+j) = x0__;")
  txt = append(txt, "\t\t\t pPrior_all(i*N__+j) = p0__;")
  txt = append(txt, "\t\t\t F_all(i*N__+j) = F__;")
  txt = append(txt, "\t\t\t A_all(i*N__+j) = A__;")
  txt = append(txt, "\t\t\t G_all(i*N__+j) = G__;")
  txt = append(txt, "\t\t }")
  txt = append(txt, "\t\t xPrior(i+1) = x0__;")
  txt = append(txt, "\t\t pPrior(i+1) = p0__;")
  txt = append(txt, "\t\t xPriorPost(2*i+2) = x0__;")
  txt = append(txt, "\t\t pPriorPost(2*i+2) = p0__;")
  
  txt = append(txt, "\n\t\t //////////// DATA-UPDATE: USING OBSERVATIONS ///////////")
  
  txt = append(txt, sprintf("\t\t y__ << %s;", paste(obs.eqs,collapse=", ")))
  txt = append(txt, "\t\t NAs__ = isNA(y__);")
  txt = append(txt, "\t\t s__ = CppAD::Integer(sum(NAs__));")
  txt = append(txt, "\t\t if( s__ > 0 ){")
  txt = append(txt, "\t\t\t ynew__ = removeNAs(s__,y__,NAs__);")
  txt = append(txt, "\t\t\t matrix<Type> E0__(s__,m__);")
  txt = append(txt, "\t\t\t E__ 	= constructE(E0__,NAs__);")
  txt = append(txt, sprintf("\t\t\t h__  = h(%s);",paste(hvars2,collapse=", ")))
  txt = append(txt, "\t\t\t e__  = ynew__ - E__*h__;")
  txt = append(txt, sprintf("\t\t\t C__  	= dhdx(%s);",paste(dhdxvars2,collapse=", ")))
  txt = append(txt,"\t\t\t C2__ 	= E__*C__;")
  txt = append(txt, sprintf("\t\t\t V__.diagonal() << %s;",paste(unlist(obsvars1),collapse=", ")))
  txt = append(txt, "\t\t\t V2__  = E__*V__*E__.transpose();")
  txt = append(txt, "\t\t\t R__ 	= C2__*p0__*C2__.transpose() + V2__;")
  txt = append(txt, "\t\t\t Ri__  = R__.inverse();")
  txt = append(txt, "\t\t\t K__ 	= p0__ * C2__.transpose() * Ri__;")
  txt = append(txt, "\t\t\t x0__  = x0__ + K__*e__;")
  txt = append(txt, "\t\t\t p0__  = (I__-K__*C2__)*p0__*(I__-K__*C2__).transpose() + K__*V2__*K__.transpose();")
  txt = append(txt, "\t\t\t nll__ += Type(0.5)*atomic::logdet(R__) + Type(0.5)*lossfunction((e__*(Ri__*e__)).sum(),tukeypars,huber_c,loss_int) + Type(0.5)*log(2*M_PI)*asDouble(s__);")
  
  txt = append(txt, "\n\t\t\t Innovation(i+1) = e__;")
  txt = append(txt, "\t\t\t InnovationCovariance(i+1) = R__;")
  txt = append(txt, "\t\t\t E_mat(i+1) = E__;")
  txt = append(txt, "\t\t\t C_mat(i+1) = C2__;")
  txt = append(txt, "\t\t\t V_mat(i+1) = V2__;")
  txt = append(txt, "\t\t\t Ri_mat(i+1) = Ri__;")
  txt = append(txt, "\t\t\t K_mat(i+1) = K__;")
  
  txt = append(txt, "\t\t }")
  txt = append(txt, "\t\t xPost(i+1) = x0__;")
  txt = append(txt, "\t\t pPost(i+1) = p0__;")
  txt = append(txt, "\t\t xPriorPost(2*i+3) = x0__;")
  txt = append(txt, "\t\t pPriorPost(2*i+3) = p0__;")
  txt = append(txt, "\t }")
  
  # Maximum-A-Posterior Estimation
  txt = append(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = append(txt, "\t //////////// MAP ESTIMATE USING PRIOR INFO ///////////")
  txt = append(txt, "\t ////////////////////////////////////////////////////////")
  txt = append(txt, "\t if(MAPbool==1){")
  txt = append(txt, "\t\t DATA_VECTOR(MAPmean);")
  txt = append(txt, "\t\t DATA_MATRIX(MAPcov);")
  txt = append(txt, sprintf("\t\t vector<Type> parvec__(%s);",length(data$pars)))
  txt = append(txt, sprintf("\t\t parvec__ << %s;",paste(names(data$pars),collapse=", ")))
  txt = append(txt, "\t\t vector<Type> pareps__ = parvec__ - MAPmean;")
  txt = append(txt, "\t\t matrix<Type> MAPcovInv__ = MAPcov.inverse();")
  txt = append(txt, "\t\t Type map_nll = Type(0.5) * atomic::logdet(MAPcov) + Type(0.5) * (pareps__ * (MAPcovInv__ * pareps__)).sum();")
  txt = append(txt, "\t\t nll__ += map_nll;")
  txt = append(txt, "\t\t REPORT(map_nll);")
  txt = append(txt, "\t\t REPORT(parvec__);")
  txt = append(txt, "\t\t REPORT(pareps__);")
  txt = append(txt, "\t }")
  
  txt = append(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = append(txt, "\t //////////// FINAL REPORTING AND RETURN //////////////")
  txt = append(txt, "\t ////////////////////////////////////////////////////////")
  txt = append(txt ,"\t REPORT(Innovation);")
  txt = append(txt ,"\t REPORT(InnovationCovariance);")
  txt = append(txt ,"\t REPORT(xPrior);")
  txt = append(txt ,"\t REPORT(xPost);")
  txt = append(txt, "\t REPORT(xPriorPost);")
  txt = append(txt, "\t REPORT(xPrior_all);")
  txt = append(txt, "\t REPORT(pPrior);")
  txt = append(txt, "\t REPORT(pPost);")
  txt = append(txt, "\t REPORT(pPriorPost);")
  txt = append(txt, "\t REPORT(pPrior_all);")
  txt = append(txt, "\t REPORT(F_all);")
  txt = append(txt, "\t REPORT(A_all);")
  txt = append(txt, "\t REPORT(G_all);")
  txt = append(txt, "\t REPORT(E_mat);")
  txt = append(txt, "\t REPORT(C_mat);")
  txt = append(txt, "\t REPORT(V_mat);")
  txt = append(txt, "\t REPORT(Ri_mat);")
  txt = append(txt, "\t REPORT(K_mat);")
  txt = append(txt, "\t return nll__;")
  txt = append(txt, "}")
  
  writeLines(txt,full_modelname)
  
  # Close file connection
  close(fileconn)
  
}
