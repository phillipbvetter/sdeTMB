write_ekf_cppfile_old = function(self, private) {

  #Initialize C++ file
  fileconn = file(paste(private$cppfile.path,".cpp",sep=""))

  txt = "#include <TMB.hpp>"
  txt = c(txt, "using namespace density;")

  ##################################################
  # CONSTRUCT AUX FUNCTIONS
  ##################################################

  txt = c(txt, "template <class Type>
vector<Type> isNA__(vector<Type> x){
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

  txt = c(txt,"template<class Type>
vector<Type> removeNAs__(int s, vector<Type> y,vector<Type> NAs){
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

  txt = c(txt, "template <class Type>
matrix<Type> constructE__(matrix<Type> E,vector<Type> p){
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

  txt = c(txt,"template<class Type>
Type lossfunction__(Type x, vector<Type> tukeypars, Type huber_c, int lossFunc){
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

  txt = c(txt, "template<class Type>
vector<Type> get_free_pars__(vector<int> mapints, int sum_mapints, vector<Type> parvec) {
	vector<Type> ans(sum_mapints);
	int j=0;
	for(int i=0;i<mapints.size();i++){
		if(mapints(i)==1){
			ans(j) = parvec(i);
			j += 1;
		}
	}
	return(ans);
}")



  ##################################################
  # CONSTRUCT DRIFT/DIFF ETC FUNCTIONS
  ##################################################

  ##################################################
  # Construct drift function

  fvars0 = sort(unique(unlist(lapply(private$diff.terms, function(x) all.vars(x$dt)))))
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars2 = fvars0
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      fvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=fvars2)
    }
  }
  for(i in 1:private$n){
    fvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0__(%s)",i-1),fvars2)
  }
  if(length(fvars2)<1){
    fvars1 = "Type a"
    fvars2 = "Type(0.0)"
  }

  temptxt = sprintf("template<class Type>\nvector<Type> f__(%s){",fvars1)
  temptxt = c(temptxt, sprintf("\tvector<Type> ans(%i);",private$n))
  for(i in 1:private$n){
    term = paste(deparse(hat2pow(private$diff.terms[[i]]$dt)),collapse="")
    temptxt = c(temptxt, sprintf("\tans(%i) = %s;",i-1,term))
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)



  ##################################################
  # Construct jacobian of drift function

  dfdxvars0 = c()
  for(i in 1:private$n){
    dfdxvars0 = c(dfdxvars0, unlist(lapply(lapply(private$state.names,function(x) Deriv::Deriv(private$diff.terms[[i]]$dt,x=x,cache.exp=FALSE)),all.vars)))
  }
  dfdxvars0 = sort(unique(dfdxvars0))
  dfdxvars1 = paste("Type", dfdxvars0, collapse=", ")
  dfdxvars2 = dfdxvars0
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      dfdxvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=dfdxvars2)
    }
  }
  for(i in 1:private$n){
    dfdxvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0__(%s)",i-1),dfdxvars2)
  }
  if(length(dfdxvars2)<1){
    dfdxvars1 = "Type a"
    dfdxvars2 = "Type(0.0)"
  }

  dfdx = matrix(0,nrow=private$n,ncol=private$n)
  for(i in 1:private$n){
    terms = lapply(private$state.names,function(x) hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt,x=x,cache.exp=FALSE)))
    for(j in 1:private$n){
      dfdx[i,j] = paste(deparse(terms[[j]]),collapse="")
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> dfdx__(%s){",dfdxvars1)
  temptxt = c(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",private$n,private$n))
  for(i in 1:private$n){
    for(j in 1:private$n){
      temptxt = c(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,dfdx[i,j]))
    }
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)



  ##################################################
  # Construct diffusion function

  gvars0 = c()
  for(i in 1:private$n){
    gvars0 = c(gvars0, unlist(lapply(private$diff.terms[[i]][-1],all.vars)))
  }
  gvars0 = sort(unique(gvars0))
  gvars1 = paste("Type", gvars0, collapse=", ")
  gvars2 = gvars0
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      gvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=gvars2)
    }
  }
  for(i in 1:private$n){
    gvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0__(%s)",i-1),gvars2)
  }

  g = matrix(0,nrow=private$n,ncol=private$ng)
  for(i in 1:private$n){
    for(j in 1:private$ng){
      term = paste(deparse(hat2pow(private$diff.terms[[i]][[j+1]])),collapse = "")
      g[i,j] = term
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> g__(%s){",gvars1)
  temptxt = c(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",private$n,private$ng))
  for(i in 1:private$n){
    for(j in 1:private$ng){
      temptxt = c(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,g[i,j]))
    }
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)



  ##################################################
  # Construct observation function

  hvars0 = unlist(lapply(lapply(private$obs.eqs.trans, function(x) x$rhs), all.vars))
  hvars0 = sort(unique(hvars0))
  hvars1 = paste("Type", hvars0, collapse=", ")
  hvars2 = hvars0
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      hvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=hvars2)
    }
  }
  for(i in 1:private$n){
    hvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0__(%s)",i-1),hvars2)
  }

  h = c()
  for(i in 1:private$m){
    term = paste(deparse(hat2pow(private$obs.eqs.trans[[i]]$rhs)),collapse = "")
    h[i] = term
  }
  temptxt = sprintf("template<class Type>\nvector<Type> h__(%s){",hvars1)
  temptxt = c(temptxt, sprintf("\tvector<Type> ans(%i);",private$m))
  for(i in 1:private$m){
    temptxt = c(temptxt, sprintf("\tans(%s) = %s;",i-1,h[i]))
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)



  ##################################################
  # Construct jacobian of observation function

  dhdxvars0 = c()
  for(i in 1:private$m){
    dhdxvars0 = c(dhdxvars0 , unlist(lapply(lapply(private$state.names, function(x) Deriv::Deriv(private$obs.eqs.trans[[i]]$rhs,x=x,cache.exp=FALSE)),all.vars)))
  }
  dhdxvars0 = sort(unique(dhdxvars0))
  dhdxvars1 = paste("Type", dhdxvars0, collapse=", ")
  dhdxvars2 = dhdxvars0
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      dhdxvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=dhdxvars2)
    }
  }
  for(i in 1:private$n){
    dhdxvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0__(%s)",i-1),dhdxvars2)
  }
  if(length(dhdxvars2)<1){
    dhdxvars1 = "Type a"
    dhdxvars2 = "Type(0.0)"
  }

  dhdx = matrix(NA,nrow=private$m,ncol=private$n)
  for(i in 1:private$m){
    terms = lapply(private$state.names, function(x) hat2pow(Deriv::Deriv(private$obs.eqs.trans[[i]]$rhs,x=x,cache.exp=FALSE)))
    for(j in 1:private$n){
      dhdx[i,j] = paste(deparse(terms[[j]]),collapse = "")
    }
  }
  temptxt = sprintf("template<class Type>\nmatrix<Type> dhdx__(%s){",dhdxvars1)
  temptxt = c(temptxt, sprintf("\tmatrix<Type> ans(%i,%i);",private$m,private$n))
  for(i in 1:private$m){
    for(j in 1:private$n){
      temptxt = c(temptxt, sprintf("\tans(%s,%s) = %s;",i-1,j-1,dhdx[i,j]))
    }
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)

  ##################################################
  # Construct observation variance diagonal function

  obsvars0 = sort(unique(unlist(lapply(private$obs.var.trans,function(x) all.vars(x$rhs)))))
  obsvars1 = paste("Type", obsvars0, collapse=", ")
  obsvars2 = obsvars0
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      obsvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=obsvars2)
    }
  }
  for(i in 1:private$n){
    obsvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0__(%s)",i-1),obsvars2)
  }
  if(length(obsvars2)<1){
    obsvars2 = "Type a"
    obsvars2 = "Type(0.0)"
  }

  temptxt = sprintf("template<class Type>\nvector<Type> obsvarFun__(%s){",obsvars1)
  temptxt = c(temptxt, sprintf("\tvector<Type> ans(%i);",private$m))
  for(i in 1:private$m){
    temptxt = c(temptxt, sprintf("\tans(%s) = %s;",i-1,paste(deparse(hat2pow(private$obs.var.trans[[i]]$rhs)),collapse="")))
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)

  ##################################################
  # BEGIN OBJECTIVE FUNCTION
  ##################################################

  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")

  ##################################################
  # Observation Vectors
  for(i in 1:length(private$obs.names)){
    txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$obs.names[i]))
  }
  # Input Vectors
  for(i in 1:length(private$input.names)){
    txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$input.names[i]))
  }
  # Initialize State
  txt = c(txt, "\t DATA_VECTOR(X0__);")
  txt = c(txt, "\t DATA_MATRIX(P0__);")
  # Time-step
  txt = c(txt, "\t DATA_VECTOR(dt__);")
  txt = c(txt, "\t DATA_IVECTOR(N__);")
  # Loss parameters
  txt = c(txt, "\t DATA_VECTOR(tukey_pars__);")
  txt = c(txt, "\t DATA_INTEGER(which_loss__);")
  txt = c(txt, "\t DATA_SCALAR(loss_c_value__);")
  # Maximum a Posterior
  txt = c(txt, "\t DATA_INTEGER(map_bool__);")

  ##################################################
  # Parameters
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("\t PARAMETER(%s);",private$parameter.names[i]))
  }

  ##################################################
  # Constants
  for (i in 1:length(private$constants)) {
    txt = c( txt , sprintf("\t DATA_SCALAR(%s);",private$constant.names[i]))
  }

  ##################################################
  # Misc
  txt = c(txt, sprintf("\n\t int n__ = %s;",private$n))
  txt = c(txt, sprintf("\t int m__ = %s;",private$m))
  txt = c(txt, "\t int s__;")
  txt = c(txt,"\t Type nll__ = 0;")

  ##################################################
  # Storage variables
  txt = c(txt, "\n\t vector<vector<Type>> xPrior(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> pPrior(t.size());")
  txt = c(txt, "\t vector<vector<Type>> xPost(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> pPost(t.size());")
  # txt = c(txt, "\t vector<vector<Type>> xPriorPost(2*t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> pPriorPost(2*t.size());")

  # txt = c(txt, "\t vector<vector<Type>> xPrior_all(k__);")
  # txt = c(txt, "\t vector<matrix<Type>> pPrior_all(k__);")

  # txt = c(txt, "\t vector<vector<Type>> F_all(t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> A_all(k__);")
  # txt = c(txt, "\t vector<matrix<Type>> G_all(k__);")

  txt = c(txt, "\t vector<vector<Type>> Innovation(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> InnovationCovariance(t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> E_mat(t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> C_mat(t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> V_mat(t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> Ri_mat(t.size());")
  # txt = c(txt, "\t vector<matrix<Type>> K_mat(t.size());")
  txt = c(txt, "\t vector<vector<Type>> H_mat(t.size());")
  txt = c(txt, "\t vector<vector<Type>> F_mat(t.size());")

  txt = c(txt, "\n\t vector<Type> x0__ = X0__;")
  txt = c(txt, "\t matrix<Type> p0__ = P0__;")
  txt = c(txt, "\t vector<Type> x1__ = x0__;")
  txt = c(txt, "\t matrix<Type> p1__ = p0__;")
  txt = c(txt, "\t xPrior(0) = x0__;")
  txt = c(txt, "\t pPrior(0) = p0__;")
  txt = c(txt, "\t xPost(0) = x0__;")
  txt = c(txt, "\t pPost(0) = p0__;")
  # txt = c(txt, "\t xPriorPost(0) = x0__;")
  # txt = c(txt, "\t xPriorPost(1) = x0__;")
  # txt = c(txt, "\t pPriorPost(0) = p0__;")
  # txt = c(txt, "\t pPriorPost(1) = p0__;")

  ##################################################
  # Initiaze variables
  txt = c(txt, "\n\t matrix<Type> C__,R__,K__,E__,C2__,V2__,Ri__,A__,G__;")
  txt = c(txt, "\t vector<Type> e0__,e__,y__(m__),NAs__,F__,ynew__,H__;")
  txt = c(txt, "\t matrix<Type> I__(n__,n__);")
  txt = c(txt, "\t I__.setIdentity();")

  ##################################################
  # Observation equation
  txt = c(txt, "\n\t matrix<Type> V__(m__,m__);")
  txt = c(txt, "\t V__.setZero();")
  obs.lhs = paste(unlist(lapply(private$obs.eqs.trans,function(x) deparse(x$form[[2]]))),"(i+1)",sep="")

  txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = c(txt, "\t //////////// MAIN FOR-LOOP ///////////")
  txt = c(txt, "\t ////////////////////////////////////////////////////////")

  ##################################################
  # Time for-loop
  txt = c(txt, "\t for(int i=0 ; i<t.size()-1 ; i++){")

  ##################################################
  # Solve Moment ODEs
  txt = c(txt, "\n\t\t //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  txt = c(txt, "\t\t for(int j=0 ; j<N__(i) ; j++){")
  txt = c(txt, sprintf("\t\t\t F__  = f__(%s);",paste(fvars2,collapse=", ")) )
  txt = c(txt, sprintf("\t\t\t A__  = dfdx__(%s);",paste(dfdxvars2,collapse=", ")) )
  txt = c(txt, sprintf("\t\t\t G__  = g__(%s);",paste(gvars2,collapse=", ")) )
  txt = c(txt,"\t\t\t x0__ = x0__ + F__ * dt__(i);")
  txt = c(txt, "\t\t\t p0__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt__(i);")
  # txt = c(txt,"\t\t\t x1__ = x0__ + F__ * dt__(i);")
  # txt = c(txt, "\t\t\t p1__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt__(i);")
  # txt = c(txt, "\t\t\t x0__ = x1__;")
  # txt = c(txt, "\t\t\t p0__ = p1__;")

  # txt = c(txt, "\n\t\t\t xPrior_all(i*N__(i)+j) = x0__;")
  # txt = c(txt, "\t\t\t pPrior_all(i*N__(i)+j) = p0__;")
  # txt = c(txt, "\t\t\t F_all(i*N__(i)+j) = F__;")
  # txt = c(txt, "\t\t\t A_all(i*N__(i)+j) = A__;")
  # txt = c(txt, "\t\t\t G_all(i*N__(i)+j) = G__;")

  txt = c(txt, "\t\t }")
  txt = c(txt, "\t\t\t F_mat(i) = F__;")
  txt = c(txt, "\t\t xPrior(i+1) = x0__;")
  txt = c(txt, "\t\t pPrior(i+1) = p0__;")
  # txt = c(txt, "\t\t xPriorPost(2*i+2) = x0__;")
  # txt = c(txt, "\t\t pPriorPost(2*i+2) = p0__;")

  ##################################################
  # Data Update
  txt = c(txt, "\n\t\t //////////// DATA-UPDATE: USING OBSERVATIONS ///////////")
  txt = c(txt, sprintf("\t\t y__ << %s;", paste(obs.lhs,collapse=", ")))
  txt = c(txt, "\t\t NAs__ = isNA__(y__);")
  txt = c(txt, "\t\t s__ = CppAD::Integer(sum(NAs__));")
  txt = c(txt, "\t\t if( s__ > 0 ){")
  txt = c(txt, "\t\t\t ynew__ = removeNAs__(s__,y__,NAs__);")
  txt = c(txt, "\t\t\t matrix<Type> E0__(s__,m__);")
  txt = c(txt, "\t\t\t E__ 	= constructE__(E0__,NAs__);")
  txt = c(txt, sprintf("\t\t\t H__  = h__(%s);",paste(hvars2,collapse=", ")))
  txt = c(txt, "\t\t\t H_mat(i+1)  = H__;")
  txt = c(txt, "\t\t\t e__  = ynew__ - E__*H__;")
  txt = c(txt, sprintf("\t\t\t C__  	= dhdx__(%s);",paste(dhdxvars2,collapse=", ")))
  txt = c(txt,"\t\t\t C2__ 	= E__*C__;")
  txt = c(txt, sprintf("\t\t\t V__.diagonal() << obsvarFun__(%s);",paste(obsvars2,collapse=" ,")))
  txt = c(txt, "\t\t\t V2__  = E__*V__*E__.transpose();")
  txt = c(txt, "\t\t\t R__ 	= C2__*p0__*C2__.transpose() + V2__;")
  txt = c(txt, "\t\t\t Ri__  = R__.inverse();")
  txt = c(txt, "\t\t\t K__ 	= p0__ * C2__.transpose() * Ri__;")
  txt = c(txt, "\t\t\t x0__  = x0__ + K__*e__;")
  txt = c(txt, "\t\t\t p0__  = (I__-K__*C2__)*p0__*(I__-K__*C2__).transpose() + K__*V2__*K__.transpose();")
  txt = c(txt, "\t\t\t nll__ += Type(0.5)*atomic::logdet(R__) + Type(0.5)*lossfunction__((e__*(Ri__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + Type(0.5)*log(2*M_PI)*asDouble(s__);")

  txt = c(txt, "\n\t\t\t Innovation(i+1) = e__;")
  txt = c(txt, "\t\t\t InnovationCovariance(i+1) = R__;")
  # txt = c(txt, "\t\t\t E_mat(i+1) = E__;")
  # txt = c(txt, "\t\t\t C_mat(i+1) = C2__;")
  # txt = c(txt, "\t\t\t V_mat(i+1) = V2__;")
  # txt = c(txt, "\t\t\t Ri_mat(i+1) = Ri__;")
  # txt = c(txt, "\t\t\t K_mat(i+1) = K__;")

  txt = c(txt, "\t\t }")
  txt = c(txt, "\t\t xPost(i+1) = x0__;")
  txt = c(txt, "\t\t pPost(i+1) = p0__;")
  # txt = c(txt, "\t\t xPriorPost(2*i+3) = x0__;")
  # txt = c(txt, "\t\t pPriorPost(2*i+3) = p0__;")
  txt = c(txt, "\t }")

  ##################################################
  # Maximum-A-Posterior
  txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = c(txt, "\t //////////// MAP ESTIMATE USING PRIOR INFO ///////////")
  txt = c(txt, "\t ////////////////////////////////////////////////////////")
  txt = c(txt, "\t if(map_bool__==1){")
  txt = c(txt, "\t\t DATA_VECTOR(map_mean__);")
  txt = c(txt, "\t\t DATA_MATRIX(map_cov__);")
  txt = c(txt, "\t\t DATA_IVECTOR(map_ints__);")
  txt = c(txt, "\t\t DATA_INTEGER(sum_map_ints__);")
  txt = c(txt, sprintf("\t\t vector<Type> parvec__(%s);",length(private$parameters)))
  txt = c(txt, "\t\t vector<Type> map_pars__;")
  txt = c(txt, sprintf("\t\t parvec__ << %s;",paste(private$parameter.names,collapse=", ")))
  txt = c(txt, sprintf("\t\t map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);"))
  txt = c(txt, "\t\t vector<Type> pars_eps__ = map_pars__ - map_mean__;")
  txt = c(txt, "\t\t matrix<Type> map_invcov__ = map_cov__.inverse();")
  txt = c(txt, "\t\t Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
  txt = c(txt, "\t\t nll__ += map_nll__;")
  txt = c(txt, "\t\t REPORT(map_nll__);")
  txt = c(txt, "\t\t REPORT(map_pars__);")
  txt = c(txt, "\t\t REPORT(pars_eps__);")
  txt = c(txt, "\t }")

  ##################################################
  # Report variables and return nll
  txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = c(txt, "\t //////////// FINAL REPORTING AND RETURN //////////////")
  txt = c(txt, "\t ////////////////////////////////////////////////////////")
  txt = c(txt ,"\t REPORT(Innovation);")
  txt = c(txt ,"\t REPORT(InnovationCovariance);")
  txt = c(txt ,"\t REPORT(xPrior);")
  txt = c(txt ,"\t REPORT(xPost);")
  txt = c(txt, "\t REPORT(pPrior);")
  txt = c(txt, "\t REPORT(pPost);")
  # txt = c(txt, "\t REPORT(xPriorPost);")
  # txt = c(txt, "\t REPORT(pPriorPost);")
  #
  # txt = c(txt, "\t REPORT(xPrior_all);")
  # txt = c(txt, "\t REPORT(pPrior_all);")

  # txt = c(txt, "\t REPORT(F_all);")
  # txt = c(txt, "\t REPORT(A_all);")
  # txt = c(txt, "\t REPORT(G_all);")
  txt = c(txt, "\t REPORT(F_mat);")

  # txt = c(txt, "\t REPORT(E_mat);")
  # txt = c(txt, "\t REPORT(C_mat);")
  # txt = c(txt, "\t REPORT(V_mat);")
  # txt = c(txt, "\t REPORT(Ri_mat);")
  txt = c(txt, "\t REPORT(H_mat);")

  txt = c(txt, "\t return nll__;")
  txt = c(txt, "}")

  writeLines(txt,fileconn)

  # Close file connection
  close(fileconn)
  #
  return(invisible(self))
}
