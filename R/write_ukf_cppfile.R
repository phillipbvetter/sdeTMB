write_ukf_cppfile = function(self, private) {

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
matrix<Type> construct_E__(matrix<Type> E,vector<Type> p){
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

  txt = c(txt,sprintf("template<class Type>
matrix<Type> construct_Xsp__(vector<Type> x0, matrix<Type> s0){
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
}",private$n,2*private$n+1)
  )

  txt = c(txt, "template<class Type>
matrix<Type> Phi__(matrix<Type> M){
  matrix<Type> K(M.col(0).size(),M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/Type(2.0);
  return K;
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

  # fvars0 = sort(unique(unlist(lapply(1:n,function(i) all.vars(diff.terms[[i]]$dt)))))
  fvars0 = sort(unique(unlist(lapply(private$diff.terms, function(x) all.vars(x$dt)))))
  fvars0_sigma = fvars0[!(fvars0 %in% private$state.names)]
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars1_sigma = paste("Type", fvars0_sigma, collapse=", ")
  fvars2 = fvars0
  fvars2_sigma = fvars0
  fvars2_sigma2 = fvars0_sigma
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      fvars2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=fvars2)
      fvars2_sigma2 = sub(pattern=sprintf("^%s$",private$input.names[i]), replacement=sprintf("%s(i)",private$input.names[i]), x=fvars2_sigma2)
    }
  }
  for(i in 1:private$n){
    fvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("__x0(%s)",i-1),fvars2)
    fvars2_sigma = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("x0(%s)",i-1),fvars2_sigma)
  }
  if(length(fvars2)<1){
    fvars1 = "Type a"
    fvars2 = "Type(0.0)"
    fvars1_sigma = "Type a"
    fvars2_sigma = "Type(0.0)"
    fvars2_sigma2 = "Type(0.0)"
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
  # Construct sigma points drift function
  ##################################################

  temptxt = sprintf("template<class Type>
matrix<Type> construct_F__(matrix<Type> Xsp, %s){
  int n = %s;
  int nn = %s;
  matrix<Type> F(n,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    F.col(i) = f__(%s);
  }
  return F;
}",fvars1_sigma,private$n,2*private$n+1,paste(fvars2_sigma,collapse=", "))
  txt = c(txt, temptxt)

  ##################################################
  # Construct diffusion function
  ##################################################

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
    gvars2 = sub(pattern=sprintf("^%s$",private$state.names[i]),replacement=sprintf("__Xsp0(0,%s)",i-1),gvars2)
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
  ##################################################

  h.allvars = unlist(lapply(lapply(private$obs.eqs.trans, function(x) x$rhs), all.vars))
  h.allvars.withoutStates = h.allvars[!(h.allvars %in% private$state.names)]
  #
  h.allvars.addedType = paste("Type", h.allvars, collapse=", ")
  h.allvars.withoutStates.addedType = paste("Type", h.allvars.withoutStates, collapse=", ")
  #
  h.allvars.addedTimeIndex.replacedStateNamesWith_x0 = h.allvars
  h.allvars.withoutStates.addedTimeIndex = h.allvars.withoutStates
  #
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      h.allvars.addedTimeIndex.replacedStateNamesWith_x0 = sub(pattern=sprintf("^%s$",private$input.names[i]),
                                                               replacement=sprintf("%s(i)",private$input.names[i]),
                                                               x=h.allvars.addedTimeIndex.replacedStateNamesWith_x0)
      h.allvars.withoutStates.addedTimeIndex = sub(pattern=sprintf("^%s$",private$input.names[i]),
                                                   replacement=sprintf("%s(i)",private$input.names[i]),
                                                   x=h.allvars.withoutStates.addedTimeIndex)
    }
  }
  for(i in 1:private$n){
    h.allvars.addedTimeIndex.replacedStateNamesWith_x0 = sub(pattern=sprintf("^%s$",private$state.names[i]),
                                                             replacement=sprintf("x0__(%s)",i-1),
                                                             x=h.allvars.addedTimeIndex.replacedStateNamesWith_x0)
  }
  if(length(h.allvars.withoutStates)<1){
    h.allvars.withoutStates.addedType = "Type a"
    h.allvars.withoutStates.addedTimeIndex = "Type(0.0)"
  }
  h.allvars.addedTimeIndex.replacedStateNamesWith_x0 = paste(h.allvars.addedTimeIndex.replacedStateNamesWith_x0, collapse=", ")
  h.allvars.withoutStates.addedTimeIndex = paste(h.allvars.withoutStates.addedTimeIndex,collapse=", ")

  h = c()
  for(i in 1:private$m){
    term = paste(deparse(hat2pow(private$obs.eqs.trans[[i]]$rhs)),collapse = "")
    h[i] = term
  }
  temptxt = sprintf("template<class Type>\nvector<Type> h__(%s){",h.allvars.addedType)
  temptxt = c(temptxt, sprintf("\tvector<Type> ans(%i);",private$m))
  for(i in 1:private$m){
    temptxt = c(temptxt, sprintf("\tans(%s) = %s;",i-1,h[i]))
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)

  ##################################################
  # Construct observation sigma points function

  temptxt = sprintf("template<class Type>
matrix<Type> construct_H__(matrix<Type> Xsp, %s){
  int nn = %s;
  int m = %s;
  matrix<Type> H(m,nn);
  vector<Type> x0__;
  for(int i=0;i<nn;i++){
    x0__ = Xsp.col(i);
    H.col(i) = h__(%s);
  }
  return H;
}",
h.allvars.withoutStates.addedType,
2*private$n+1,
private$m,
h.allvars.addedTimeIndex.replacedStateNamesWith_x0
  )
  txt = c(txt, temptxt)

  ##################################################
  # Construct observation variance diagonal function

  hvar.allvars = sort(unique(unlist(lapply(private$obs.var.trans,function(x) all.vars(x$rhs)))))
  hvar.allvars.withoutStates = hvar.allvars[!(hvar.allvars %in% private$state.names)]
  #
  hvar.allvars.addedType = paste("Type", hvar.allvars, collapse=", ")
  hvar.allvars.withoutStates.addedType = paste("Type", hvar.allvars.withoutStates, collapse=", ")
  #
  hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0 = hvar.allvars
  hvar.allvars.withoutStates.addedTimeIndex = hvar.allvars.withoutStates
  #
  if(length(private$input.names)>0){
    for(i in 1:length(private$input.names)){
      hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0 = sub(pattern=sprintf("^%s$",private$input.names[i]),
                                                                  replacement=sprintf("%s(i)",private$input.names[i]),
                                                                  x=hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0)
      hvar.allvars.withoutStates.addedTimeIndex = sub(pattern=sprintf("^%s$",private$input.names[i]),
                                                      replacement=sprintf("%s(i)",private$input.names[i]),
                                                      x=hvar.allvars.withoutStates.addedTimeIndex)
    }
  }
  for(i in 1:private$n){
    hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0 = sub(pattern=sprintf("^%s$",private$state.names[i]),
                                                                replacement=sprintf("x0__(%s)",i-1),
                                                                x=hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0)
  }
  if(length(hvar.allvars.withoutStates)<1){
    hvar.allvars.withoutStates.addedType = "Type a"
    hvar.allvars.withoutStates.addedTimeIndex = "Type(0.0)"
  }
  hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0 = paste(hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0, collapse=", ")
  hvar.allvars.withoutStates.addedTimeIndex = paste(hvar.allvars.withoutStates.addedTimeIndex,collapse=", ")

  # TO DO
  # When calling we give Xsp.col(0) and then rename states to x0(0),x0(1) etc...
  temptxt = sprintf("template<class Type>\nvector<Type> obsvar_diagonal__(%s){",hvar.allvars.addedType)
  temptxt = c(temptxt, sprintf("\tvector<Type> ans(%i);",private$m))
  for(i in 1:private$m){
    temptxt = c(temptxt, sprintf("\tans(%s) = %s;",i-1,paste(deparse(hat2pow(private$obs.var.trans[[i]]$rhs)),collapse="")))
  }
  temptxt = c(temptxt, "\treturn ans;\n}")
  txt = c(txt, temptxt)

  temptxt = sprintf("template<class Type>
vector<Type> obsvar_diagonal_usingXsp__(matrix<Type> Xsp, %s){
	vector<Type> ans;
  vector<Type> x0__ = Xsp.col(0);
	ans = obsvar_diagonal__(%s);
	return ans;
}",
hvar.allvars.withoutStates.addedType,
hvar.allvars.addedTimeIndex.replacedStateNamesWith_x0
  )
  txt = c(txt, temptxt)

  ##################################################
  # Write TMB objective function
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
  txt = c(txt, sprintf("\t int nn__ = %s;",2*private$n+1))
  txt = c(txt, sprintf("\t int m__ = %s;",private$m))
  txt = c(txt, "\t int s__;")
  txt = c(txt, "\t Type nll__ = 0;")

  ##################################################
  # Storage variables
  txt = c(txt, "\n\t vector<vector<Type>> xPrior(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> pPrior(t.size());")
  txt = c(txt, "\t vector<vector<Type>> xPost(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> pPost(t.size());")
  txt = c(txt, "\t vector<vector<Type>> Innovation(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> InnovationCovariance(t.size());")

  txt = c(txt, "\t xPrior(0) = X0__;")
  txt = c(txt, "\t pPrior(0) = P0__;")
  txt = c(txt, "\t xPost(0) = X0__;")
  txt = c(txt, "\t pPost(0) = P0__;")

  ##################################################
  # Initiaze variables
  txt = c(txt, "\n\t vector<Type> y__(m__);")
  txt = c(txt, "\t vector<Type> NAs__,ynew__,e__,X1__;")
  txt = c(txt, "\t matrix<Type> Xsp0__,S0__,G__, F__, S0Inv__, M__, S0PhiM__, Frhs1__(n__,nn__), Frhs0__, Frhs__;")
  txt = c(txt, "\t matrix<Type> E__, H__, Syy__, SyyInv__, Sxy__, K__, P1__;")


  ##################################################
  # Observation equation
  txt = c(txt,"\t matrix<Type> V__(m__,m__);")
  txt = c(txt,"\t V__.setZero();")
  obs.lhs = paste(unlist(lapply(private$obs.eqs.trans,function(x) deparse(x$form[[2]]))),"(i+1)",sep="")

  ##################################################
  # Weights
  txt = c(txt, "\n\t Type lambda__ = 3 - n__, weights__;")
  txt = c(txt, "\t vector<Type> Wm__(nn__);")
  txt = c(txt, "\t matrix<Type> Wmm__(nn__,nn__), Wm_diag__(nn__,nn__), I__(nn__,nn__), W__;")
  txt = c(txt, "\t weights__ = Type(1.0)/(Type(2.0)*(lambda__ + n__));")
  txt = c(txt, "\t Wm__.fill(weights__);")
  txt = c(txt, "\t Wm__(0) = lambda__/(lambda__ + n__);")
  txt = c(txt, "\t for(int i=0; i<nn__ ; i++){")
  txt = c(txt, "\t\t Wmm__.col(i) = Wm__;")
  txt = c(txt, "\t }")
  txt = c(txt, "\t Wm_diag__.setZero();")
  txt = c(txt, "\t Wm_diag__.diagonal() = Wm__;")
  txt = c(txt, "\t I__.setIdentity();")
  txt = c(txt, "\t W__ = (I__ - Wmm__) * Wm_diag__ * (I__ - Wmm__).transpose();\n")

  txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
  txt = c(txt, "\t //////////// MAIN FOR-LOOP ///////////")
  txt = c(txt, "\t ////////////////////////////////////////////////////////")

  txt = c(txt, "\n\t S0__ = P0__.llt().matrixL();")
  txt = c(txt, "\t Xsp0__ = construct_Xsp__(X0__,S0__);")

  ##################################################
  # Time for-loop
  txt = c(txt, "\n\t for(int i=0 ; i<t.size()-1 ; i++){")

  ##################################################
  # Solve Moment ODEs
  txt = c(txt, "\n\t\t //////////// Time-Update ///////////")
  txt = c(txt, "\t\t for(int j=0 ; j<N__(i) ; j++){")
  txt = c(txt, sprintf("\t\t\t F__  = construct_F__(Xsp0__,%s);",paste(fvars2_sigma2,collapse=", ")))
  txt = c(txt, sprintf("\t\t\t G__  = g__(%s);",paste(gvars2,collapse=", ")) )
  txt = c(txt, "\t\t\t S0Inv__ = S0__.inverse();")
  txt = c(txt, "\t\t\t M__ = S0Inv__ * (Xsp0__ * W__ * F__.transpose() + F__ * W__ * Xsp0__.transpose() + G__*G__.transpose()) * S0Inv__.transpose();")
  txt = c(txt, "\t\t\t S0PhiM__ = S0__ * Phi__(M__);")
  txt = c(txt, "\t\t\t Frhs1__.block(0,1,n__,n__) = S0PhiM__;")
  txt = c(txt, "\t\t\t Frhs1__.block(0,n__+1,n__,n__) = -S0PhiM__;")
  txt = c(txt, "\t\t\t Frhs0__ = (F__*Wm__).replicate(1,nn__);")
  txt = c(txt, "\t\t\t Frhs__ = Frhs0__ + sqrt(3.0) * Frhs1__;")
  txt = c(txt, "\t\t\t Xsp0__ += Frhs__ * dt__(i);")
  txt = c(txt, "\t\t\t S0__ = ((Xsp0__ - Xsp0__.col(0).replicate(1,nn__))/sqrt(Type(3.0))).block(0,1,n__,n__);")
  txt = c(txt, "\t\t\t};")
  txt = c(txt, "\t\t P1__ = S0__ * S0__.transpose();")
  txt = c(txt, "\t\t xPrior(i+1) = Xsp0__.col(0);;")
  txt = c(txt, "\t\t pPrior(i+1) = P1__;")

  txt = c(txt, "\n\t\t //////////// Time-Update ///////////")
  txt = c(txt, sprintf("\t\t y__ << %s;", paste(obs.lhs,collapse=", ")))
  txt = c(txt, "\t\t NAs__ = isNA__(y__);")
  txt = c(txt, "\t\t s__ = CppAD::Integer(sum(NAs__));")
  txt = c(txt, "\t\t if( s__ > 0 ){")
  txt = c(txt, "\t\t\t ynew__ = removeNAs__(s__,y__,NAs__);")
  txt = c(txt, "\t\t\t matrix<Type> E0__(s__,m__);")
  txt = c(txt, "\t\t\t E__ 	= construct_E__(E0__,NAs__);")
  txt = c(txt, sprintf("\t\t\t H__ = construct_H__(Xsp0__,%s);",h.allvars.withoutStates.addedTimeIndex))
  txt = c(txt, "\t\t\t e__  = ynew__ - E__ * (H__ * Wm__);")
  txt = c(txt, sprintf("\t\t\t V__.diagonal() << obsvar_diagonal_usingXsp__(Xsp0__, %s);",hvar.allvars.withoutStates.addedTimeIndex))
  txt = c(txt, "\t\t\t Syy__  = E__ * (H__ * W__ * H__.transpose() + V__) * E__.transpose();")
  txt = c(txt, "\t\t\t Sxy__  = Xsp0__ * W__ * H__.transpose() * E__.transpose();")
  txt = c(txt, "\t\t\t SyyInv__  = Syy__.inverse();")
  txt = c(txt, "\t\t\t K__ = Sxy__ * SyyInv__;")
  txt = c(txt, "\t\t\t X1__ = Xsp0__ * Wm__ + K__ * e__;")
  txt = c(txt, "\t\t\t P1__ = S0__ * S0__.transpose() - K__ * Syy__ * K__.transpose();")
  txt = c(txt, "\t\t\t S0__ = P1__.llt().matrixL();")
  txt = c(txt, "\t\t\t Xsp0__ = construct_Xsp__(X1__,S0__);")
  # txt = c(txt, "\t\t\t nll__ += 0.5*atomic::logdet(Syy__) + 0.5*(e__*(SyyInv__*e__)).sum() + 0.5*log(2*M_PI)*asDouble(s__);")
  txt = c(txt, "\t\t\t nll__ += Type(0.5)*atomic::logdet(Syy__) + 0.5*lossfunction__((e__*(SyyInv__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + Type(0.5)*log(2*M_PI)*asDouble(s__);")
  txt = c(txt, "\t\t\t Innovation(i+1) = e__;")
  txt = c(txt, "\t\t\t InnovationCovariance(i+1) = Syy__;")
  txt = c(txt, "\t\t };")
  txt = c(txt, "\t\t xPost(i+1) = Xsp0__.col(0);")
  txt = c(txt, "\t\t pPost(i+1) = P1__;")
  txt = c(txt, "\t };")

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
  txt = c(txt ,"\t REPORT(xPrior);")
  txt = c(txt, "\t REPORT(pPrior);")
  txt = c(txt,"\t REPORT(xPost);")
  txt = c(txt,"\t REPORT(pPost);")
  txt = c(txt,"\t REPORT(Innovation);")
  txt = c(txt,"\t REPORT(InnovationCovariance);")

  txt = c(txt, "\t return nll__;")
  txt = c(txt, "};")

  writeLines(txt,fileconn)

  # Close file connection
  close(fileconn)

  return(invisible(self))
}
