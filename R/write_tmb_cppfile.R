write_tmb_cppfile = function(self, private) {

  #Initialize C++ file
  fileconn = file(paste(private$cppfile.path,".cpp",sep=""))

  txt = "#include <TMB.hpp>"
  txt = c(txt, "using namespace density;")

  ##################################################
  # CONSTRUCT AUX FUNCTIONS
  ##################################################

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

  writeLines(txt,fileconn)

  # CREATE TIMEDEPINPUTS FROM INPUTS AND STATES (STATES ARE RANDOM EFFECTS
  # PARAMETERS IN THE TMB FRAMEWORK)

  timedepInput = c(private$input.names,private$state.names)

  ##################################################
  # CONSTRUCT DRIFT/DIFF ETC FUNCTIONS
  ##################################################

  ##################################################
  # Construct drift function

  fvars0 = sort(unique(unlist(lapply(private$diff.terms, function(x) all.vars(x$dt)))))
  fvars1 = paste("Type", fvars0, collapse=", ")
  fvars2 = fvars0
  input.names = private$input.names
  for(i in 1:length(private$input.names)){
    fvars2 = sub(pattern=sprintf("^%s$",input.names[i]), replacement=sprintf("%s(i)",input.names[i]), x=fvars2)
  }
  state.names = private$state.names
  for(i in 1:length(private$state.names)){
    fvars2 = sub(pattern=sprintf("^%s$",state.names[i]), replacement=sprintf("%s(Nc__(i)+j)",state.names[i]), x=fvars2)
  }
  temptxt = sprintf("template<class Type>\nvector<Type> f__(%s){",fvars1)
  temptxt = c(temptxt, sprintf("\tvector<Type> f(%i);",private$n))
  for(i in 1:private$n){
    term = paste(deparse(hat2pow(private$diff.terms[[i]]$dt)),collapse="")
    temptxt = c(temptxt, sprintf("\tf(%i) = %s;",i-1,term))
  }
  temptxt = c(temptxt, "\treturn(f);\n}")
  txt = c(txt, temptxt)

  ##################################################
  # Construct diffusion matrix function

  gvars0 = c()
  for(i in 1:private$n){
    gvars0 = c(gvars0, unlist(lapply(private$diff.terms[[i]][-1],all.vars)))
  }
  gvars0 = sort(unique(gvars0))
  gvars1 = paste("Type", gvars0, collapse=", ")
  gvars2 = gvars0
  input.names = private$input.names
  for(i in 1:length(private$input.names)){
    gvars2 = sub(pattern=sprintf("^%s$",input.names[i]), replacement=sprintf("%s(i)",input.names[i]), x=gvars2)
  }
  state.names = private$state.names
  for(i in 1:length(private$state.names)){
    gvars2 = sub(pattern=sprintf("^%s$",state.names[i]), replacement=sprintf("%s(Nc__(i)+j)",state.names[i]), x=gvars2)
  }

  temptxt = sprintf("template<class Type>\nmatrix<Type> g__(%s){",gvars1)
  temptxt = c(temptxt, sprintf("\tmatrix<Type> G(%i,%i);",private$n,private$ng))
  for(i in 1:private$n){
    for(j in 1:private$ng){
      term = paste(deparse(hat2pow(private$diff.terms[[i]][[j+1]])),collapse="")
      temptxt = c( temptxt , sprintf("\tG(%i,%i) = %s;",i-1,j-1,term))
    }
  }
  temptxt = c(temptxt, "\treturn(G);\n}")
  txt = c(txt, temptxt)


  ##################################################
  # Construct observation stuff

  hvars0 = sort(unique(unlist(lapply(lapply(private$obs.eqs.trans, function(x) x$rhs), all.vars))))
  hvars1 = paste("Type", hvars0, collapse=", ")
  hvars2 = hvars0
  input.names = private$input.names
  for(i in 1:length(private$input.names)){
    hvars2 = sub(pattern=sprintf("^%s$",input.names[i]), replacement=sprintf("%s(j)",input.names[i]), x=hvars2)
  }
  state.names = private$state.names
  for(i in 1:length(private$state.names)){
    hvars2 = sub(pattern=sprintf("^%s$",state.names[i]), replacement=sprintf("%s(Nc__(j))",state.names[i]), x=hvars2)
  }


  ##################################################
  # Construct observation variance diagonal function

  obsvars0 = sort(unique(unlist(lapply(private$obs.var.trans,function(x) all.vars(x$rhs)))))
  obsvars1 = paste("Type", obsvars0, collapse=", ")
  obsvars2 = obsvars0
  input.names = private$input.names
  for(i in 1:length(private$input.names)){
    obsvars2 = sub(pattern=sprintf("^%s$",input.names[i]), replacement=sprintf("%s(i)",input.names[i]), x=obsvars2)
  }
  state.names = private$state.names
  for(i in 1:length(private$state.names)){
    obsvars2 = sub(pattern=sprintf("^%s$",state.names[i]), replacement=sprintf("%s(Nc__(i))",state.names[i]), x=obsvars2)
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

  # Initialize objective function
  txt = c(txt, "\n//////////// objective function ///////////")
  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()
  {")

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

  # State Vectors
  txt = c(txt, "\n//// state random effect vectors ////")
  for(i in 1:length(private$state.names)){
    txt = c(txt, sprintf("\t PARAMETER_VECTOR(%s);",private$state.names[i]))
  }

  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "\t DATA_VECTOR(X0__);")
  txt = c(txt, "\t DATA_MATRIX(P0__);")

  # Time-step
  txt = c(txt, "\n//// time-step ////")
  txt = c(txt, "\t DATA_VECTOR(dt__);")
  txt = c(txt, "\t DATA_IVECTOR(N__);")
  txt = c(txt, "\t DATA_IVECTOR(Nc__);")

  # Iobs
  txt = c(txt, "\n//// iobs vectors ////")
  for (i in 1:private$m) {
    nam = paste("iobs_",private$obs.names[i],sep="")
    txt = c(txt, sprintf("\t DATA_IVECTOR(%s);",nam))
  }

  # Loss parameters
  # txt = c(txt, "\t DATA_VECTOR(tukey_pars__);")
  # txt = c(txt, "\t DATA_INTEGER(which_loss__);")
  # txt = c(txt, "\t DATA_SCALAR(loss_c_value__);")

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
  txt = c( txt , "\t DATA_INTEGER(ng__);")

  # Likelihood
  # txt = c(txt,sprintf("\n\t int n__ = %s;",private$n))
  # txt = c(txt,sprintf("\t int ng__ = %s;",private$ng))
  # txt = c(txt,sprintf("\t int m__ = %s;",private$m))
  # txt = c(txt,"\t Type nll__ = 0;")

  # Initiaze variables
  txt = c(txt, "\n\t //////////// initialize variables ///////////")
  txt = c(txt,"\t Type nll__ = 0;")
  txt = c(txt, "\t vector<Type> F__(n__);")
  txt = c(txt, "\t vector<Type> Xi__(n__);")
  txt = c(txt, "\t vector<Type> Xip1__(n__);")
  txt = c(txt, "\t vector<Type> Z__(n__);")
  txt = c(txt, "\t matrix<Type> G__(n__,ng__);")
  txt = c(txt, "\t matrix<Type> V__(n__,n__);")
  txt = c(txt, "\t matrix<Type> I__(n__,n__);")
  txt = c(txt, "\t I__.setIdentity();")
  txt = c(txt, "\t I__ *= 1e-8;")

  # Forward simulation and likelihood contribution from states
  txt = c(txt, "\n\t //////////// MAIN LOOP OVER TIME POINTS ///////////")

  txt = c(txt, "\t for(int i=0 ; i<t.size()-1 ; i++){")
  #
  txt = c(txt, "\t\t for(int j=0 ; j<N__(i) ; j++){")
  #
  txt = c(txt, sprintf("\t\t\t F__ = f__(%s);",paste(fvars2,collapse=", ")))
  txt = c(txt, sprintf("\t\t\t G__ = g__(%s);",paste(gvars2,collapse=", ")))
  txt = c(txt, sprintf("\t\t\t Xi__ << %s;",paste(private$state.names,"(Nc__(i)+j)",collapse=", ",sep="")))
  txt = c(txt, sprintf("\t\t\t Xip1__ << %s;",paste(private$state.names,"(Nc__(i)+j+1)",collapse=", ",sep="")))
  txt = c(txt, sprintf("\t\t\t Z__ = Xip1__ - Xi__ - F__ * dt__(i);"))
  txt = c(txt, sprintf("\t\t\t V__ = (G__ * G__.transpose() + I__) * dt__(i);"))
  txt = c(txt, "\t\t\t nll__ += MVNORM(V__)(Z__);")
  #
  txt = c(txt, "\t\t }")
  #
  txt = c(txt, "\t }")

  # Datta Update
  txt = c(txt, "\n\t\t //////////// DATA-UPDATE ///////////")
  #
  txt = c(txt, "\t matrix<Type> varDiag__(m__,t.size());")
  txt = c(txt, "\t for(int i=0 ; i<t.size() ; i++){")
  txt = c(txt, sprintf("\t\t varDiag__.col(i) = obsvarFun__(%s);",paste(obsvars2,collapse=", ")))
  txt = c(txt, "\t }")
  #
  for(i in 1:private$m){
    nam = paste("iobs_",private$obs.names[i],sep="")
    txt = c(txt, sprintf("\t for(int i=0 ; i<%s.size() ; i++){",nam))
    txt = c(txt, sprintf("\t\t int j = %s(i);",nam))
    txt = c(txt, sprintf("\t\t nll__ -= dnorm(%s,%s,sqrt(varDiag__.col(j)(%s)),true);",paste(private$obs.names[i],"(j)",sep=""),
                         hvars2[i],i-1))
    txt = c(txt, "\t }")
  }

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
  txt = c(txt, "\t return nll__;")
  txt = c(txt, "}")

  # Write cpp function and close file connection
  writeLines(txt,fileconn)
  close(fileconn)

  # Return
  return(invisible(self))
}
