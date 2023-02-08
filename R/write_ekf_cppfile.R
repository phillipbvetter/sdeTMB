write_ekf_cppfile = function(self, private) {

  #Initialize C++ file
  fileconn = file(paste(private$cppfile.path,".cpp",sep=""))

  txt = "#include <TMB.hpp>"
  txt = c(txt, "using namespace density;")

  ##################################################
  # CONSTRUCT AUX FUNCTIONS
  ##################################################

  # Function for NA-checking
  txt = c(txt, "\n//////////// helper fun: find NA locations in vector ///////////")
  txt = c(txt, "template <class Type>
vector<Type> is_not_na(vector<Type> x){
	vector<Type> y(x.size());
	y.fill(Type(1.0));
    for(int i=0; i<x.size(); i++){
		if( R_IsNA(asDouble(x(i))) ){
			y(i) = Type(0.0);
		}
	}
	return y;
}")

  # Implements Tukey and Huber loss functions
  txt = c(txt, "\n//////////// helper fun: extract non-NAs from vector ///////////")
  txt = c(txt,"template<class Type>
vector<Type> remove_nas(vector<Type> data_vector, int number_of_datas, vector<Type> na_bool){
  int ii = 0;
	vector<Type> y_reduced(number_of_datas);
	for(int i=0; i < data_vector.size(); i++){
		if(na_bool(i) == Type(1.0)){
			y_reduced(ii) = data_vector(i);
			ii++;
		}
	}
  return y_reduced;
}")

  # Implements Tukey and Huber loss functions
  txt = c(txt, "\n//////////// helper fun: construct permutation matrix ///////////")
  txt = c(txt, "template <class Type>
matrix<Type> construct_permutation_matrix(int number_of_datas, int m, vector<Type> na_bool){
	matrix<Type> E(number_of_datas,m);
	E.setZero();
	/**/
	int j=0;
	for(int i=0; i < m; i++){
		if(na_bool(i) == Type(1.0)){ /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
			E(j,i) = Type(1.0);
			j += 1;
		}
	}
	return E;
}")

  # Implements Tukey and Huber loss functions
  txt = c(txt, "\n//////////// Loss function ///////////")
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

  # Helper function for MAP estimation
  txt = c(txt, "\n//////////// Map estimation helper ///////////")
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

#   # Allows import of list of integer matrices (for permutation matrices in data update)
#   txt = c(txt, "\n//////////// list of matrices ///////////")
#   txt = c(txt, "template<class Type>
# struct LOIM_t : vector<matrix<Type>> {
#   LOIM_t(SEXP x){  /* x = list of vectors passed from R */
#     (*this).resize(LENGTH(x));
#     for(int i=0; i<LENGTH(x); i++){
#       SEXP sm = VECTOR_ELT(x, i);
#       (*this)(i) = asMatrix<Type>(sm);
#     }
#   }
# };")
#
#   # Allows import of list of numeric vectors
#   txt = c(txt, "\n//////////// list of vectors ///////////")
#   txt = c(txt, "template<class Type>
# struct LONV_t : vector<vector<Type>> {
#   LONV_t(SEXP x){  /* x = list of vectors passed from R */
#     (*this).resize(LENGTH(x));
#     for(int i=0; i<LENGTH(x); i++){
#       SEXP sm = VECTOR_ELT(x, i);
#       (*this)(i) = asVector<Type>(sm);
#     }
#   }
# };")



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

  # temptxt = sprintf("template<class Type>\nvector<Type> obsvarFun__(%s){",obsvars1)
  # temptxt = c(temptxt, sprintf("\tvector<Type> ans(%i);",private$m))
  # for(i in 1:private$m){
  #   temptxt = c(temptxt, sprintf("\tans(%s) = %s;",i-1,paste(deparse(hat2pow(private$obs.var.trans[[i]]$rhs)),collapse="")))
  # }
  # temptxt = c(temptxt, "\treturn ans;\n}")
  # txt = c(txt, temptxt)
  temptxt = sprintf("template<class Type>\nmatrix<Type> obsvarFun__(%s){",obsvars1)
  temptxt = c(temptxt, sprintf("\tmatrix<Type> V(%i,%i);",private$m,private$m))
  for(i in 1:private$m){
    temptxt = c(temptxt, sprintf("\tV(%s,%s) = %s;",i-1,i-1,paste(deparse(hat2pow(private$obs.var.trans[[i]]$rhs)),collapse="")))
  }
  temptxt = c(temptxt, "\treturn V;\n}")
  txt = c(txt, temptxt)

  ##################################################
  # BEGIN OBJECTIVE FUNCTION
  ##################################################

  # Initialize objective function
  txt = c(txt, "\n//////////// objective function ///////////")
  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")

  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  for(i in 1:length(private$obs.names)){
    txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$obs.names[i]))
  }

  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  for(i in 1:length(private$input.names)){
    txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$input.names[i]))
  }

  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "\t DATA_VECTOR(X0__);")
  txt = c(txt, "\t DATA_MATRIX(P0__);")

  # Time-step
  txt = c(txt, "\t DATA_VECTOR(dt__);")
  txt = c(txt, "\t DATA_IVECTOR(N__);")

  # Loss parameters
  txt = c(txt, "\n//// loss parameters ////")
  txt = c(txt, "\t DATA_VECTOR(tukey_pars__);")
  txt = c(txt, "\t DATA_INTEGER(which_loss__);")
  txt = c(txt, "\t DATA_SCALAR(loss_c_value__);")

  # Maximum a Posterior
  txt = c(txt, "\n//// map estimation ////")
  txt = c(txt, "\t DATA_INTEGER(map_bool__);")

  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("\t PARAMETER(%s);",private$parameter.names[i]))
  }

  # Constants
  txt = c(txt, "\n//// constants ////")
  for (i in 1:length(private$constants)) {
    txt = c( txt , sprintf("\t DATA_SCALAR(%s);",private$constant.names[i]))
  }

  # system size
  txt = c(txt, "\n//// system size ////")
  txt = c( txt , "\t DATA_INTEGER(n__);")
  txt = c( txt , "\t DATA_INTEGER(m__);")

  # Storage variables
  txt = c(txt, "\n//////////// storage variables ///////////")
  txt = c(txt, "\t vector<vector<Type>> xPrior(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> pPrior(t.size());")
  txt = c(txt, "\t vector<vector<Type>> xPost(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> pPost(t.size());")
  txt = c(txt, "\t vector<vector<Type>> Innovation(t.size());")
  txt = c(txt, "\t vector<matrix<Type>> InnovationCovariance(t.size());")

  txt = c(txt, "\n//////////// set initial value ///////////")
  txt = c(txt, "\t vector<Type> x0__ = X0__;")
  txt = c(txt, "\t matrix<Type> p0__ = P0__;")
  txt = c(txt, "\t xPrior(0) = X0__;")
  txt = c(txt, "\t xPost(0) = X0__;")
  txt = c(txt, "\t pPrior(0) = P0__;")
  txt = c(txt, "\t pPost(0) = P0__;")

  # Initiaze variables
  txt = c(txt, "\n\t //////////// initialize variables ///////////")
  txt = c(txt, "\t int s__;")
  txt = c(txt, "\t Type half_log2PI = Type(0.5)*log(2*M_PI);")
  txt = c(txt, "\t Type nll__ = 0;")
  txt = c(txt, "\t vector<Type> data_vector__(m__),na_bool__,e__,y__,F__,H__;")
  txt = c(txt, "\t matrix<Type> C__,R__,K__,E__,V__,Ri__,A__,G__;")

  # Identity Matrix
  txt = c(txt, "\n\t //////////// identity matrix ///////////")
  txt = c(txt, "\t matrix<Type> I__(n__,n__);")
  txt = c(txt, "\t I__.setIdentity();")

  # Main For-Loop
  txt = c(txt, "\n\t //////////// MAIN LOOP OVER TIME POINTS ///////////")

  # Time for-loop
  txt = c(txt, "\t for(int i=0 ; i<t.size()-1 ; i++){")

  # Solve Moment ODEs
  txt = c(txt, "\n\t\t //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  txt = c(txt, "\t\t for(int j=0 ; j<N__(i) ; j++){")
  #
  txt = c(txt, sprintf("\t\t\t F__  = f__(%s);",paste(fvars2,collapse=", ")) )
  txt = c(txt, sprintf("\t\t\t A__  = dfdx__(%s);",paste(dfdxvars2,collapse=", ")) )
  txt = c(txt, sprintf("\t\t\t G__  = g__(%s);",paste(gvars2,collapse=", ")) )
  #
  txt = c(txt,"\t\t\t x0__ = x0__ + F__ * dt__(i);")
  txt = c(txt, "\t\t\t p0__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt__(i);")
  #
  txt = c(txt, "\t\t }")
  #
  txt = c(txt, "\t\t xPrior(i+1) = x0__;")
  txt = c(txt, "\t\t pPrior(i+1) = p0__;")

  # Data Update
  txt = c(txt, "\n\t\t //////////// DATA-UPDATE ///////////")
  obs.lhs = paste(unlist(lapply(private$obs.eqs.trans,function(x) deparse(x$form[[2]]))),"(i+1)",sep="")
  txt = c(txt, sprintf("\t\t data_vector__ << %s;", paste(obs.lhs,collapse=", ")))
  txt = c(txt, "\t\t na_bool__ = is_not_na(data_vector__);")
  txt = c(txt, "\t\t s__ = CppAD::Integer(sum(na_bool__));")
  # Start if statement for observations
  txt = c(txt, "\t\t if( s__ > 0 ){")
  txt = c(txt, "\t\t\t y__  = remove_nas(data_vector__, s__, na_bool__);")
  txt = c(txt, "\t\t\t E__  = construct_permutation_matrix(s__, m__, na_bool__);")
  txt = c(txt, sprintf("\t\t\t H__  = h__(%s);",paste(hvars2,collapse=", ")))
  txt = c(txt, sprintf("\t\t\t C__  = E__ * dhdx__(%s);", paste(dhdxvars2,collapse=", ")))
  txt = c(txt, "\t\t\t e__  = y__ - E__ * H__;")
  txt = c(txt, sprintf("\t\t\t V__  = E__ * obsvarFun__(%s) * E__.transpose();", paste(obsvars2,collapse=" ,")))
  txt = c(txt, "\t\t\t R__  = C__ * p0__ * C__.transpose() + V__;")
  txt = c(txt, "\t\t\t Ri__ = R__.inverse();")
  txt = c(txt, "\t\t\t K__ 	= p0__ * C__.transpose() * Ri__;")
  txt = c(txt, "\t\t\t x0__ = x0__ + K__*e__;")
  txt = c(txt, "\t\t\t p0__ = (I__ - K__ * C__) * p0__ * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();")
  txt = c(txt, "\t\t\t nll__ += Type(0.5)*atomic::logdet(R__) + Type(0.5)*lossfunction__((e__*(Ri__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + half_log2PI * asDouble(s__);")
  txt = c(txt, "\t\t\t Innovation(i+1) = e__;")
  txt = c(txt, "\t\t\t InnovationCovariance(i+1) = R__;")
  txt = c(txt, "\t\t }")
  # end if statement for observations
  txt = c(txt, "\t\t xPost(i+1) = x0__;")
  txt = c(txt, "\t\t pPost(i+1) = p0__;")
  #
  txt = c(txt, "\t }")


  # Maximum-A-Posterior
  txt = c(txt, "\n\t\t //////////// MAP CONTRIBUTION ///////////")
  #
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


  # Report variables and return nll
  txt = c(txt, "\n\t //////////// Return/Report //////////////")
  #
  txt = c(txt ,"\t REPORT(Innovation);")
  txt = c(txt ,"\t REPORT(InnovationCovariance);")
  txt = c(txt ,"\t REPORT(xPrior);")
  txt = c(txt ,"\t REPORT(xPost);")
  txt = c(txt, "\t REPORT(pPrior);")
  txt = c(txt, "\t REPORT(pPost);")
  txt = c(txt, "\t return nll__;")
  txt = c(txt, "}")

  # Write cpp function and close file connection
  writeLines(txt,fileconn)
  close(fileconn)

  # Return
  return(invisible(self))
}
