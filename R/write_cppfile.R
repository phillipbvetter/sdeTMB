
########## GENERAL HELPER FUNCTIONS ###############
############################################################
write_helper_cppfunctions = function(){
  
  txt = c()
  
  # Find function finds indices of NAs in a vector
  newtxt = "\n//////////// FIND NA INDICES IN VECTOR ///////////
  template <class Type>
  vector<Type> is_not_na(vector<Type> x){
    vector<Type> y(x.size());
    y.fill(Type(1.0));
      for(int i=0; i<x.size(); i++){
        if( R_IsNA(asDouble(x(i))) ){
          y(i) = Type(0.0);
        }
      }
    return y;
  }"
  txt = c(txt,newtxt)
  
  # This function removes NAs from a vector
  newtxt = "\n//////////// REMOVE NA'S FROM VECTOR ///////////
  template<class Type>
  vector<Type> remove_nas__(vector<Type> obsVec, int number_of_nonNA_observations, vector<Type> is_not_na_vector){
    int ii = 0;
    vector<Type> y_reduced(number_of_nonNA_observations);
      for(int i=0; i < obsVec.size(); i++){
        if(is_not_na_vector(i) == Type(1.0)){
          y_reduced(ii) = obsVec(i);
          ii++;
        }
      }
    return y_reduced;
  }"
  txt = c(txt,newtxt)
  
  # Construct Permutation Matrix E for Kalman Filter NA-removal
  newtxt = "\n//////////// helper fun: construct permutation matrix ///////////
  template <class Type>
  matrix<Type> construct_permutation_matrix(int number_of_nonNA_observations, int number_of_obs_eqs, vector<Type> is_not_na_vector){
	  matrix<Type> E(number_of_nonNA_observations, number_of_obs_eqs);
	  E.setZero();
	  int j=0;
	  for(int i=0; i < number_of_obs_eqs; i++){
      /*if p(i) is 1 then include by setting 1 in diagonal of matrix*/
		  if(is_not_na_vector(i) == Type(1.0)){
        E(j,i) = Type(1.0);
			  j += 1;
      }
	  }
	  return E;
  }"
  txt = c(txt,newtxt)
  
  # Implements Tukey and Huber loss functions
  newtxt = "\n//////////// loss function ///////////
  template<class Type>
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
    return loss;
  }"
  txt = c(txt,newtxt)
  
  # Helper function for MAP estimation
  newtxt = "\n//////////// MAP estimation helper ///////////
  template<class Type>
  vector<Type> get_free_pars__(vector<int> mapints, int sum_mapints, vector<Type> parvec) {
	  vector<Type> ans(sum_mapints);
	  int j=0;
	  for(int i=0;i<mapints.size();i++){
		  if(mapints(i)==1){
			  ans(j) = parvec(i);
			  j += 1;
		  }
	  }
	  return ans;
  }"
  txt = c(txt, newtxt)
  
  return(txt)
}

########## UNSCENTED KALMAN FILTER FUNCTIONS ###############
############################################################
write_ukf_helper_functions = function(self, private){
  
  txt = c()
  
  # Construct sigma points function
  newtxt = "\n//////////// UKF sigma points ///////////
  template<class Type>
  matrix<Type> construct_Xsp__(vector<Type> x0, matrix<Type> s0){
    int n = %s;
    int nn = %s;
    matrix<Type> Xsp(n, nn);
    vector<Type> Si;
    Xsp.col(0) = x0;
    for(int i=1; i<n+1; i++){
      Si = s0.col(i-1);
      Xsp.col(i) = x0 + sqrt(1+n) * Si;
      Xsp.col(i+n) = x0 - sqrt(1+n) * Si;
    }
    return Xsp;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, 2*private$number.of.states+1)
  txt = c(txt, newtxt)
  
  # Construct Phi Matrix for Sigma Points
  newtxt = "\n//////////// UKF function ///////////
  template<class Type>
  matrix<Type> Phi__(matrix<Type> M){
    matrix<Type> K(M.col(0).size(), M.row(0).size());
    K.setZero();
    K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
    K.diagonal() = K.diagonal()/Type(2.0);
    return K;
  }"
  txt = c(txt, newtxt)
  
  # Construct sigma points drift function
  newtxt = "\n//////////// UKF sigma points drift function ///////////
  template<class Type>
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
  }"
  newtxt = sprintf(newtxt, fvars1_sigma, private$n, 2*private$n+1, paste(fvars2_sigma,collapse=", "))
  txt = c(txt, newtxt)
  
  # return
  return(txt) 
}

write_ekf_functions = function(self, private){
  
  txt = c()
  
  # Create substitution translation list
  obsList = lapply(seq_along(private$obs.names), function(id) substitute(obsVec(i),list(i=as.numeric(id-1))))
  names(obsList) = private$obs.names
  parList = lapply(seq_along(private$parameter.names), function(id) substitute(parVec(i),list(i=as.numeric(id-1))))
  names(parList) = private$parameter.names
  stateList = lapply(seq_along(private$state.names), function(id) substitute(stateVec(i),list(i=as.numeric(id-1))))
  names(stateList) = private$state.names
  inputList = lapply(seq_along(private$input.names), function(id) substitute(inputVec(i),list(i=as.numeric(id-1))))
  names(inputList) = private$input.names
  subsList = c(obsList, parList, stateList, inputList)
  
  ##################################################
  # drift
  ##################################################
  
  # Perform substitution of parameters, inputs and states
  f = sapply( seq_along(private$state.names),
              function(i){
                drift.term = hat2pow(private$diff.terms[[i]]$dt)
                new.drift.term = do.call(substitute, list(drift.term, subsList))
                sprintf("f__(%i) = %s;",i-1,deparse1(new.drift.term))
              })
  
  newtxt = "\n//////////// drift function //////////
  template<class Type>
  vector<Type> f__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> f__(%s);
    %s
    return f__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, paste(f,collapse="\n\t\t"))
  txt = c(txt,newtxt)
  
  ##################################################
  # drift jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  jac.f = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(Deriv::Deriv(private$diff.terms[[i]]$dt,x=private$state.names[j], cache.exp=F))
      new.term = do.call(substitute, list(term, subsList))
      jac.f = c(jac.f, sprintf("dfdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// jacobian of drift function ///////////
  template<class Type>
  matrix<Type> dfdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dfdx__(%s, %s);
    %s
    return dfdx__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.states, paste(jac.f, collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # diffusion
  ##################################################
  
  # calculate all the terms and substitute variables
  g = c()
  for(i in seq_along(private$state.names)){
    for(j in seq_along(private$diff.processes[-1])){
      term = hat2pow(private$diff.terms[[i]][[j+1]])
      new.term = do.call(substitute, list(term, subsList))
      g = c(g, sprintf("g__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  newtxt = "\n//////////// diffusion function ///////////
  template<class Type>
  matrix<Type> g__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> g__(%s, %s);
    %s
    return g__;
  }"
  newtxt = sprintf(newtxt, private$number.of.states, private$number.of.diffusions, paste(g,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation
  ##################################################
  
  # calculate all the terms and substitute variables
  h = sapply(seq_along(private$obs.names), 
             function(i){
               term = hat2pow(private$obs.eqs.trans[[i]]$rhs)
               new.term = do.call(substitute, list(term, subsList))
               sprintf("h__(%s) = %s;",i-1, deparse1(new.term))
             }) 
  
  newtxt = "\n//////////// observation function ///////////
  template<class Type>
  vector<Type> h__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> h__(%s);
    %s
    return h__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(h,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation jacobian
  ##################################################
  
  # calculate all the terms and substitute variables
  jac.h = c()
  for(i in seq_along(private$obs.names)){
    for(j in seq_along(private$state.names)){
      term = hat2pow(private$diff.terms.obs[[i]][[j]])
      new.term = do.call(substitute, list(term, subsList))
      jac.h = c(jac.h, sprintf("dhdx__(%s, %s) = %s;",i-1, j-1, deparse1(new.term)))
    }
  }
  
  newtxt = "\n//////////// jacobian of obs function ///////////
  template<class Type>
  matrix<Type> dhdx__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> dhdx__(%s, %s);
    %s
    return dhdx__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, private$number.of.states, paste(jac.h,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # observation variance
  ##################################################
  
  obs.var = lapply(seq_along(private$obs.var.trans), 
                   function(i) {
                     term = hat2pow(private$obs.var.trans[[i]]$rhs)
                     new.term = do.call(substitute, list(term, subsList))
                     sprintf("hvar__(%s) = %s;", i-1, deparse1(new.term))
                   })
  newtxt = "\n//////////// observation variance matrix function ///////////
  template<class Type>
  vector<Type> hvar__(vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    vector<Type> hvar__(%s);
    %s
    return hvar__;
  }"
  newtxt = sprintf(newtxt, private$number.of.observations, paste(obs.var,collapse="\n\t\t"))
  txt = c(txt, newtxt)
  
  ##################################################
  # 1-Step MOMENT DIFFERENTIAL EQUATIONS
  ##################################################
  
  newtxt = "\n//////////// 1-step f moment ODE ///////////
  template<class Type>
  matrix<Type> cov_ode_1step(matrix<Type> covMat, vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec){
    matrix<Type> cov_ode_1step = dfdx__(stateVec, parVec, inputVec) * covMat + covMat * dfdx__(stateVec, parVec, inputVec).transpose() + g__(stateVec, parVec, inputVec) * g__(stateVec, parVec, inputVec).transpose();
    return cov_ode_1step;
  }"
  txt = c(txt, newtxt)
  
  ##################################################
  # ODE SOLVER
  ##################################################
  
  newtxt = "\n//////////// ODE SOLVER ///////////
  template<class Type>
  struct ode_integration {
  
	  vector<Type> X1;
	  matrix<Type> P1;
  
	  ode_integration(matrix<Type> covMat, vector<Type> stateVec, vector<Type> parVec, vector<Type> inputVec, vector<Type> dinputVec, Type dt, int ode_solver){
      /*Initial State and Cov Values*/
      vector<Type> X0 = stateVec;
      matrix<Type> P0 = covMat;
  
		  /*Forward Euler*/
		  if(ode_solver == 1){
  
		   X1 = X0 + f__(stateVec, parVec, inputVec) * dt;
       P1 = P0 + cov_ode_1step(covMat, stateVec, parVec, inputVec) * dt;
  
		  /*4th Order Runge-Kutta 4th*/
		  } else if (ode_solver == 2){
  
		   vector<Type> k1, k2, k3, k4;
		   matrix<Type> c1, c2, c3, c4;
  
		   /*1. Approx Slope at Initial Point*/
       k1 = f__(stateVec, parVec, inputVec);
       c1 = cov_ode_1step(covMat, stateVec, parVec, inputVec);
  
		   /*2. First Approx Slope at Midpoint*/
       // inputVec += 0.5 * dinputVec;
       stateVec = X0 + 0.5 * dt * k1;
       covMat   = P0 + 0.5 * dt * c1;
       k2       = f__(stateVec, parVec, inputVec); 
       c2       = cov_ode_1step(covMat, stateVec, parVec, inputVec);        
  
		   /*3. Second Approx Slope at Midpoint*/
       stateVec = X0 + 0.5 * dt * k2;
       covMat   = P0 + 0.5 * dt * c2;
       k3       = f__(stateVec, parVec, inputVec); 
       c3       = cov_ode_1step(covMat, stateVec, parVec, inputVec);
  
		   /*4. Approx Slope at End Point*/
       // inputVec += 0.5 * dinputVec;
       stateVec = X0 + dt * k3;
       covMat   = P0 + dt * c3;
       k4       = f__(stateVec, parVec, inputVec); 
       c4       = cov_ode_1step(covMat, stateVec, parVec, inputVec);
  
		   /*ODE UPDATE*/
		   X1 = X0 + dt * (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
		   P1 = P0 + dt * (c1 + 2.0*c2 + 2.0*c3 + c4)/6.0;
  
		 } else {
		 /*nothing*/
		 }
		}
  };"
  txt = c(txt, newtxt)
  
  # return
  return(txt)
}

write_ekf_estimate = function(self, private){
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  txt = c(txt, "DATA_MATRIX(obsMat)")

  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  txt = c(txt, "DATA_MATRIX(inputMat)")
  
  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "DATA_VECTOR(stateVec);")
  txt = c(txt, "DATA_MATRIX(covMat);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  
  # Loss parameters
  txt = c(txt, "\n//// loss parameters ////")
  txt = c(txt, "DATA_VECTOR(tukey_loss_parameters);")
  txt = c(txt, "DATA_INTEGER(loss_function);")
  txt = c(txt, "DATA_SCALAR(loss_threshold_value);")
  
  # Maximum a Posterior
  txt = c(txt, "\n//// map estimation ////")
  txt = c(txt, "DATA_INTEGER(MAP_bool);")
  
  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("PARAMETER(%s);",private$parameter.names[i]))
  }
  
  # system size
  txt = c(txt, "\n//// system size ////")
  txt = c( txt , "DATA_INTEGER(number_of_state_eqs);")
  txt = c( txt , "DATA_INTEGER(number_of_obs_eqs);")
  txt = c( txt , "int tsize = inputMat.col(0).size();")
  
  # state, par, input, obs
  txt = c(txt, "\n//// state, par, input, obs vectors ////")
  txt = c(txt, sprintf("vector<Type> inputVec(%s), dinputVec(%s), obsVec(%s), parVec(%s);", 
                       private$number.of.inputs, 
                       private$number.of.inputs,
                       private$number.of.observations,
                       private$number.of.pars
                       ))
  txt = c(txt, sprintf("parVec << %s;", paste(private$parameter.names,collapse=", ")))
  
  # Storage variables
  txt = c(txt, "\n//////////// storage variables ///////////")
  txt = c(txt, "vector<vector<Type>> xPrior(tsize), xPost(tsize), Innovation(tsize);")
  txt = c(txt, "vector<matrix<Type>> pPrior(tsize), pPost(tsize), InnovationCovariance(tsize);")
  
  txt = c(txt, "\n//////////// set initial value ///////////")
  txt = c(txt, "xPrior(0) = stateVec, xPost(0) = stateVec;")
  txt = c(txt, "pPrior(0) = covMat, pPost(0) = covMat;")
  
  # Initiaze variables
  txt = c(txt, "\n //////////// initialize variables ///////////")
  txt = c(txt, "int number_of_nonNA_observations;")
  txt = c(txt, "Type half_log2PI = Type(0.5) * log(2*M_PI);")
  txt = c(txt, "vector<Type> is_not_na_obsVec, e__, y__, F__, H__;")
  txt = c(txt, "matrix<Type> C__, R__, K__, E__, V__, Ri__, A__, G__;")
  
  # Identity Matrix
  txt = c(txt, "\n//////////// identity matrix ///////////")
  txt = c(txt, "matrix<Type> I__(number_of_state_eqs, number_of_state_eqs), V0__(number_of_obs_eqs, number_of_obs_eqs);")
  txt = c(txt, "I__.setIdentity();")
  
  ########## <MAIN LOOP>  ##########
  txt = c(txt, "\n //////////// MAIN LOOP OVER TIME POINTS ///////////")
  
  # Time for-loop
  txt = c(txt, "for(int i=0 ; i < tsize - 1 ; i++){")
  # Set inputs
  txt = c(txt, "inputVec = inputMat.row(i);")
  txt = c(txt, "dinputVec = (inputMat.row(i+1) - inputMat.row(i))/ode_timesteps(i);")
  # Solve Moment ODEs
  txt = c(txt, "\n //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  txt = c(txt, "for(int j=0 ; j < ode_timesteps(i) ; j++){")
  txt = c(txt, "ode_integration<Type> odelist = {covMat, stateVec, parVec, inputVec, dinputVec, ode_timestep_size(i), ode_solver};")
  txt = c(txt, "stateVec = odelist.X1;")
  txt = c(txt, "covMat = odelist.P1;")
  # txt = c(txt, "inputVec += dinputVec;")
  txt = c(txt, "}")
  txt = c(txt, "xPrior(i+1) = stateVec;")
  txt = c(txt, "pPrior(i+1) = covMat;")
  
  # Data Update
  txt = c(txt, "\n //////////// DATA-UPDATE ///////////")
  
  # Set observation indices (i+1)
  txt = c(txt, "obsVec = obsMat.row(i+1);")
  
  # Check the number of NAs in the observation vector
  txt = c(txt, "is_not_na_obsVec = is_not_na(obsVec);")
  txt = c(txt, "number_of_nonNA_observations = CppAD::Integer(sum(is_not_na_obsVec));")
  
  ########## <OBS IF STATEMENT>  ##########
  txt = c(txt, "if( number_of_nonNA_observations > 0 ){")
  txt = c(txt, "inputVec = inputMat.row(i+1);")
  txt = c(txt, "y__ = remove_nas__(obsVec, number_of_nonNA_observations, is_not_na_obsVec);")
  txt = c(txt, "E__ = construct_permutation_matrix(number_of_nonNA_observations, number_of_obs_eqs, is_not_na_obsVec);")
  txt = c(txt, "H__ = h__(stateVec, parVec, inputVec);")
  txt = c(txt, "C__ = E__ * dhdx__(stateVec, parVec, inputVec);")
  txt = c(txt, "e__ = y__ - E__ * H__;")
  txt = c(txt, "V0__.diagonal() << hvar__(stateVec, parVec, inputVec);")
  txt = c(txt, "V__ = E__ * V0__ * E__.transpose();")
  txt = c(txt, "R__ = C__ * covMat * C__.transpose() + V__;")
  txt = c(txt, "Ri__ = R__.inverse();")
  txt = c(txt, "K__ = covMat * C__.transpose() * Ri__;")
  txt = c(txt, "stateVec = stateVec + K__*e__;")
  txt = c(txt, "covMat = (I__ - K__ * C__) * covMat * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();")
  txt = c(txt, "nll__ += Type(0.5) * atomic::logdet(R__) + Type(0.5) * lossfunction__((e__*(Ri__*e__)).sum(), tukey_loss_parameters, loss_threshold_value, loss_function) + half_log2PI * asDouble(number_of_nonNA_observations);")
  txt = c(txt, "Innovation(i+1) = e__;")
  txt = c(txt, "InnovationCovariance(i+1) = R__;")
  txt = c(txt, "}")
  ########## </OBS IF STATEMENT>  ##########
  
  # Store posterior values
  txt = c(txt, "xPost(i+1) = stateVec;")
  txt = c(txt, "pPost(i+1) = covMat;")
  
  txt = c(txt, "}")
  ########## </MAIN LOOP>  ##########
  
  # Maximum-A-Posterior
  txt = c(txt, "\n//////////// MAP CONTRIBUTION ///////////")
  txt = c(txt, "if(MAP_bool == 1){")
  txt = c(txt, "DATA_VECTOR(map_mean__);")
  txt = c(txt, "DATA_MATRIX(map_cov__);")
  txt = c(txt, "DATA_IVECTOR(map_ints__);")
  txt = c(txt, "DATA_INTEGER(sum_map_ints__);")
  txt = c(txt, "vector<Type> map_pars__;")
  txt = c(txt, "map_pars__ = get_free_pars__(map_ints__, sum_map_ints__, parVec);")
  txt = c(txt, "vector<Type> pars_eps__ = map_pars__ - map_mean__;")
  txt = c(txt, "matrix<Type> map_invcov__ = map_cov__.inverse();")
  txt = c(txt, "Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
  txt = c(txt, "nll__ += map_nll__;")
  txt = c(txt, "REPORT(map_nll__);")
  txt = c(txt, "REPORT(map_pars__);")
  txt = c(txt, "REPORT(pars_eps__);")
  txt = c(txt, "}")
  
  # Report variables
  txt = c(txt, "\n//////////// Report //////////////")
  txt = c(txt, "REPORT(Innovation);")
  txt = c(txt, "REPORT(InnovationCovariance);")
  txt = c(txt, "REPORT(xPrior);")
  txt = c(txt, "REPORT(xPost);")
  txt = c(txt, "REPORT(pPrior);")
  txt = c(txt, "REPORT(pPost);")
  
  # return
  return(txt)
}

write_ekf_predict = function(self, private){
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  for(i in 1:length(private$obs.names)){
    txt = c(txt, sprintf("DATA_VECTOR(%s);",private$obs.names[i]))
  }
  
  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  txt = c(txt, "DATA_MATRIX(inputMat)")
  # for(i in 1:length(private$input.names)){
  # txt = c(txt, sprintf("DATA_VECTOR(%s);",private$input.names[i]))
  # }
  
  # Initialize State
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "DATA_VECTOR(stateVec);")
  txt = c(txt, "DATA_MATRIX(covMat);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  
  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in 1:length(private$parameters)){
    txt = c(txt, sprintf("PARAMETER(%s);",private$parameter.names[i]))
  }
  
  # system size
  txt = c(txt, "\n//// system size ////")
  txt = c( txt , "DATA_INTEGER(number_of_state_eqs);")
  txt = c( txt , "DATA_INTEGER(number_of_obs_eqs);")
  txt = c(txt, "int number_of_state_eqs_squared = number_of_state_eqs * number_of_state_eqs;")
  
  # state, par, input, obs
  txt = c(txt, "\n//// state, par, input, obs vectors ////")
  txt = c(txt, sprintf("vector<Type> parVec(%s);", private$number.of.pars))
  txt = c(txt, sprintf("parVec << %s;", paste(private$parameter.names,collapse=", ")))
  txt = c(txt, sprintf("vector<Type> obsVec(%s);", private$number.of.observations))
  txt = c(txt, sprintf("vector<Type> inputVec(%s), dinputVec(%s);", 
                       private$number.of.inputs, 
                       private$number.of.inputs))
  
  # Prediction variables
  txt = c( txt , "DATA_INTEGER(last_pred_index);")
  txt = c( txt , "DATA_INTEGER(k_step_ahead);")
  
  # Prediction storage
  # index i storage
  txt = c(txt, "vector<matrix<Type>> xk__(last_pred_index);")
  txt = c(txt, "vector<matrix<Type>> pk__(last_pred_index);")
  txt = c(txt, "xk__.fill(matrix<Type>(k_step_ahead+1, number_of_state_eqs));")
  txt = c(txt, "pk__.fill(matrix<Type>(k_step_ahead+1, number_of_state_eqs_squared));")
  # temp index k storage
  txt = c(txt, "matrix<Type> xk_temp__(k_step_ahead+1, number_of_state_eqs);")
  txt = c(txt, "array<Type> pk_temp__(number_of_state_eqs, number_of_state_eqs, k_step_ahead+1);")
  
  # in each k-step-ahead im saving k state vectors, and k covariance matrices temporarily.
  # these can be saved into a matrix and stored in a vector of index i...
  
  # Initiaze variables
  txt = c(txt, "\n //////////// initialize variables ///////////")
  txt = c(txt, "int number_of_nonNA_observations;")
  txt = c(txt, "vector<Type> is_not_na_obsVec, e__, y__, F__, H__;")
  txt = c(txt, "matrix<Type> C__, R__, K__, E__, V__, Ri__, A__, G__;")
  
  # Identity Matrix
  txt = c(txt, "\n//////////// identity matrix ///////////")
  txt = c(txt, "matrix<Type> I__(number_of_state_eqs, number_of_state_eqs), V0__(number_of_obs_eqs, number_of_obs_eqs);")
  txt = c(txt, "I__.setIdentity();")
  
  # # Main For-Loop
  txt = c(txt, "\n //////////// MAIN LOOP OVER TIME POINTS ///////////")
  txt = c(txt, "for(int i=0 ; i < last_pred_index ; i++){")
  
  txt = c(txt,"xk_temp__.row(0) = stateVec;")
  txt = c(txt,"pk_temp__.col(0) = covMat;")
  
  # K-step-ahead for-loop
  txt = c(txt, "\n //////////// K-STEP-AHEAD LOOP ///////////")
  txt = c(txt, "for(int k=0 ; k < k_step_ahead ; k++){")
  
  # Set input indices (i+k)
  txt = c(txt, "inputVec = inputMat.row(i+k);")
  txt = c(txt, "dinputVec = (inputMat.row(i+k+1) - inputMat.row(i+k))/ode_timesteps(i+k);")
  
  ########################################################
  # TIME UPDATE
  ########################################################
  txt = c(txt, "\n //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////")
  txt = c(txt, "for(int j=0 ; j < ode_timesteps(i+k) ; j++){")
  txt = c(txt, "ode_integration<Type> odelist = {covMat, stateVec, parVec, inputVec, ode_timestep_size(i+k), ode_solver};")
  txt = c(txt, "stateVec = odelist.X1;")
  txt = c(txt, "covMat = odelist.P1;")
  # txt = c(txt, sprintf("inputVec(0) += ode_timestep_size(i+k);"))
  txt = c(txt, "inputVec += dinputVec;") #update first order time inputs
  txt = c(txt, "}")
  # Save each prediction from k = {0, 1, ..., k_step_ahead-1, k_step_ahead}
  txt = c(txt, "xk_temp__.row(k+1) = stateVec;")
  txt = c(txt, "pk_temp__.col(k+1) = covMat;")
  txt = c(txt, "}")
  
  # Save all these k-step predictions as 
  txt = c(txt, "xk__(i) = xk_temp__;")
  txt = c(txt, "for(int kk=0 ; kk < k_step_ahead+1 ; kk++){")
  txt = c(txt, "pk__(i).row(kk) = vector<Type>(pk_temp__.col(kk).transpose());")
  txt = c(txt, "}")
  
  # Extract 1-step predictions for next i-loop
  txt = c(txt,"stateVec = xk_temp__.row(1);")
  txt = c(txt,"covMat = pk_temp__.col(1).matrix();")
  
  ########################################################
  # DATA UPDATE
  ########################################################
  txt = c(txt, "\n //////////// DATA-UPDATE ///////////")
  
  # Set observation indices (i+1)
  obs.ids = paste0( private$obs.names, "(i+1)", collapse=", ")
  txt = c(txt, sprintf("obsVec << %s;", obs.ids))
  
  # Check the number of NAs in the observation vector
  txt = c(txt, "is_not_na_obsVec = is_not_na(obsVec);")
  txt = c(txt, "number_of_nonNA_observations = CppAD::Integer(sum(is_not_na_obsVec));")
  
  ########## <OBSERVATION IF STATEMENT>  ##########
  txt = c(txt, "if( number_of_nonNA_observations > 0 ){")
  # input.ids = paste0( private$input.names, "(i+1)", collapse=", ") # Indice input indices (i+1)
  # txt = c(txt, sprintf("inputVec << %s;", input.ids))
  txt = c(txt, sprintf("inputVec = inputMat.row(i+1);"))
  txt = c(txt, "y__ = remove_nas__(obsVec, number_of_nonNA_observations, is_not_na_obsVec);") # Remove NAs in the observation vector
  txt = c(txt, "E__ = construct_permutation_matrix(number_of_nonNA_observations, number_of_obs_eqs, is_not_na_obsVec);") # Permutation matrix to filter out NA observation indices
  txt = c(txt, "H__ = h__(stateVec, parVec, inputVec);")
  txt = c(txt, "C__ = E__ * dhdx__(stateVec, parVec, inputVec);")
  txt = c(txt, "e__ = y__ - E__ * H__;")
  txt = c(txt, "V0__.diagonal() << hvar__(stateVec, parVec, inputVec);")
  txt = c(txt, "V__ = E__ * V0__ * E__.transpose();")
  txt = c(txt, "R__ = C__ * covMat * C__.transpose() + V__;")
  txt = c(txt, "Ri__ = R__.inverse();")
  txt = c(txt, "K__ = covMat * C__.transpose() * Ri__;")
  txt = c(txt, "stateVec = stateVec + K__*e__;")
  txt = c(txt, "covMat = (I__ - K__ * C__) * covMat * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();")
  txt = c(txt, "}")
  # # Report variables
  txt = c(txt, "\n //////////// Return and Report //////////////")
  txt = c(txt ,"REPORT(xk__);")
  txt = c(txt ,"REPORT(pk__);")
  txt = c(txt, "}")
  
  # return
  return(txt)
}

write_laplace_estimate = function(self, private){
  
  txt = c()
  
  # Observation Vectors
  txt = c(txt, "\n//// observations ////")
  for(i in seq_along(private$obs.names)){
    txt = c(txt, sprintf("DATA_VECTOR(%s);",private$obs.names[i]))
  }
  
  # Input Vectors
  txt = c(txt, "\n//// inputs ////")
  for(i in seq_along(private$input.names)){
    txt = c(txt, sprintf("DATA_VECTOR(%s);",private$input.names[i]))
  }
  
  # Random Effects State Vectors
  txt = c(txt, "\n//// state random effect vectors ////")
  for(i in seq_along(private$state.names)){
    txt = c(txt, sprintf("PARAMETER_VECTOR(%s);",private$state.names[i]))
  }
  
  # Initialize State
  # Thee are used, but can be used to evaluate initial likelihood at some point?
  txt = c(txt, "\n//// initial state ////")
  txt = c(txt, "DATA_VECTOR(stateVec);")
  txt = c(txt, "DATA_MATRIX(covMat);")
  
  # Time-step
  txt = c(txt, "DATA_VECTOR(ode_timestep_size);")
  txt = c(txt, "DATA_IVECTOR(ode_timesteps);")
  txt = c(txt, "DATA_IVECTOR(ode_cumsum_timesteps);")
  
  # Iobs
  txt = c(txt, "\n//// iobs vectors ////")
  for (i in seq_along(private$obs.names)) {
    nam = paste("iobs_",private$obs.names[i],sep="")
    txt = c(txt, sprintf("DATA_IVECTOR(%s);",nam))
  }
  
  # Maximum a Posterior
  txt = c(txt, "\n//// map estimation ////")
  txt = c(txt, "DATA_INTEGER(MAP_bool);")
  
  # Parameters
  txt = c(txt, "\n//// parameters ////")
  for(i in seq_along(private$parameters)){
    txt = c(txt, sprintf("PARAMETER(%s);",private$parameter.names[i]))
  }
  
  # system size
  txt = c(txt, "\n//// system size ////")
  txt = c(txt, "DATA_INTEGER(number_of_state_eqs);")
  txt = c(txt, "DATA_INTEGER(number_of_obs_eqs);")
  txt = c(txt, "DATA_INTEGER(number_of_diffusions);")
  
  # state, par, input, obs
  txt = c(txt, "\n//// state, par, input, obs vectors ////")
  txt = c(txt, sprintf("vector<Type> parVec(%s);", private$number.of.pars))
  txt = c(txt, sprintf("parVec << %s;", paste(private$parameter.names,collapse=", ")))
  txt = c(txt, sprintf("vector<Type> inputVec(%s);", private$number.of.inputs))
  
  # Initiaze variables
  txt = c(txt, "\n //////////// initialize variables ///////////")
  txt = c(txt, "vector<Type> F__(number_of_state_eqs);")
  txt = c(txt, "vector<Type> stateVec0(number_of_state_eqs);")
  txt = c(txt, "vector<Type> stateVec1(number_of_state_eqs);")
  txt = c(txt, "vector<Type> Z__(number_of_state_eqs);")
  txt = c(txt, "matrix<Type> G__(number_of_state_eqs, number_of_diffusions);")
  txt = c(txt, "matrix<Type> V__(number_of_state_eqs, number_of_state_eqs);")
  txt = c(txt, "matrix<Type> I__(number_of_state_eqs, number_of_state_eqs);")
  txt = c(txt, "matrix<Type> varDiag__(number_of_obs_eqs, t.size());")
  txt = c(txt, "matrix<Type> varFun__(number_of_obs_eqs, t.size());")
  txt = c(txt, "I__.setIdentity();")
  txt = c(txt, "I__ *= 1e-8;")
  
  # Forward simulation and likelihood contribution from states
  txt = c(txt, "\n\t //////////// MAIN LOOP OVER TIME POINTS ///////////")
  
  txt = c(txt, "\t for(int i=0 ; i < t.size()-1 ; i++){")
  txt = c(txt, "\t\t for(int j=0 ; j < ode_timesteps(i) ; j++){")
  txt = c(txt, sprintf("stateVec0 << %s;", paste0(private$state.names,"(ode_cumsum_timesteps(i)+j)",collapse=", ")))
  txt = c(txt, sprintf("stateVec1 << %s;", paste0(private$state.names,"(ode_cumsum_timesteps(i)+j+1)",collapse=", ")))
  txt = c(txt, sprintf("inputVec << %s;", paste0(private$input.names,"(ode_cumsum_timesteps(i)+j)",collapse=", ")))
  txt = c(txt, "F__ = f__(stateVec0, parVec, inputVec);")
  txt = c(txt, "G__ = g__(stateVec0, parVec, inputVec);")
  txt = c(txt, "Z__ = stateVec1 - stateVec0 - F__ * ode_timestep_size(i);")
  txt = c(txt, "V__ = (G__ * G__.transpose() + I__) * ode_timestep_size(i);")
  txt = c(txt, "nll__ += MVNORM(V__)(Z__);")
  txt = c(txt, "}")
  txt = c(txt, "}")
  
  # Data Update
  txt = c(txt, "\n\t\t //////////// DATA-UPDATE ///////////")
  # We calculate the observation equation and variance for all observations (also NAs)
  txt = c(txt, "for(int i=0 ; i < t.size() ; i++){")
  txt = c(txt, sprintf("stateVec0 << %s;", paste0(private$state.names,"(i)",collapse=", ")))
  txt = c(txt, sprintf("inputVec << %s;", paste0(private$input.names,"(i)",collapse=", ")))
  txt = c(txt, "varDiag__.col(i) = hvar__(stateVec0, parVec, inputVec);")
  txt = c(txt, "varFun__.col(i) = h__(stateVec0, parVec, inputVec);")
  txt = c(txt, "}")
  
  # This is a marginal for-loop for each observation
  for(i in seq_along(private$obs.names)){
    nam = paste0("iobs_",private$obs.names[i])
    txt = c(txt, sprintf("for(int i=0 ; i < %s.size() ; i++){",nam))
    txt = c(txt, sprintf("int j = %s(i);",nam))
    txt = c(txt, sprintf("nll__ -= dnorm(%s, varFun__.col(j)(%s), sqrt(varDiag__.col(j)(%s)), true);",
                         paste0(private$obs.names[i],"(j)"),
                         i-1,
                         i-1))
    txt = c(txt, "}")
  }
  
  # Maximum-A-Posterior
  txt = c(txt, "\n//////////// MAP CONTRIBUTION ///////////")
  txt = c(txt, "if(MAP_bool == 1){")
  txt = c(txt, "DATA_VECTOR(map_mean__);")
  txt = c(txt, "DATA_MATRIX(map_cov__);")
  txt = c(txt, "DATA_IVECTOR(map_ints__);")
  txt = c(txt, "DATA_INTEGER(sum_map_ints__);")
  txt = c(txt, "vector<Type> map_pars__;")
  txt = c(txt, sprintf("map_pars__ = get_free_pars__(map_ints__, sum_map_ints__, parVec);"))
  txt = c(txt, "vector<Type> pars_eps__ = map_pars__ - map_mean__;")
  txt = c(txt, "matrix<Type> map_invcov__ = map_cov__.inverse();")
  txt = c(txt, "Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
  txt = c(txt, "nll__ += map_nll__;")
  txt = c(txt, "REPORT(map_nll__);")
  txt = c(txt, "REPORT(map_pars__);")
  txt = c(txt, "REPORT(pars_eps__);")
  txt = c(txt, "}")
  
  # return
  return(txt)
}

write_cppfile = function(self, private) {
  
  #Initialize C++ file
  fileconn = file(paste0(private$cppfile.path,".cpp"))
  
  # Include TMB header
  txt = "#include <TMB.hpp>"
  
  # density namespace for special functions
  newtxt = "using namespace density;"
  txt = c(txt, newtxt)
  
  ##################################################
  # GENERAL HELPER FUNCTIONS
  ##################################################
  
  newtxt = write_helper_cppfunctions()
  txt = c(txt, newtxt)
  
  ##################################################
  # HELPER FUNCTIONS FOR UNSCENTED KALMAN FILTER
  ##################################################
  
  # newtxt = write_ukf_helper_functions(self,private)
  # txt = c(txt,newtxt)
  
  ####################################################################
  # FUNCTIONS FOR DRIFT, DIFFUSION, OBSERVATION, OBS VARIANCE, ETC...
  ####################################################################
  
  newtxt = write_ekf_functions(self, private)
  txt = c(txt, newtxt)
  
  ##################################################
  # BEGIN OBJECTIVE FUNCTION
  ##################################################
  
  # Initialize objective function
  txt = c(txt, "\n//////////// objective function ///////////")
  txt = c(txt,"template<class Type>\nType objective_function<Type>::operator() ()\n{")
  
  txt = c(txt, "DATA_INTEGER(estimation_method);")
  txt = c(txt, "DATA_INTEGER(prediction_bool);")
  txt = c(txt, "DATA_INTEGER(ode_solver);")
  txt = c(txt, "Type nll__ = 0;")
  
  ##################################################
  # EKF OBJECTIVE FUNCTION
  ##################################################
  
  txt = c(txt, "\n//////////// EKF METHOD ///////////")
  txt = c(txt, "if (estimation_method == 1 && prediction_bool == 0) {")
  
  newtxt = write_ekf_estimate(self, private)
  txt = c(txt, newtxt)
  
  txt = c(txt, "} else if (estimation_method == 1 && prediction_bool == 1) {")
  
  # newtxt = write_ekf_predict(self, private)
  # txt = c(txt, newtxt)
  
  txt = c(txt, "} else if (estimation_method == 2 && prediction_bool == 0) {")
  
  # newtxt = write_ukf_estimate(self, private)
  # txt = c(txt, newtxt)
  
  txt = c(txt, "} else if (estimation_method == 2 && prediction_bool == 1) {")
  
  # newtxt = write_ukf_predict(self, private)
  # txt = c(txt, newtxt)
  
  txt = c(txt, "} else if (estimation_method == 3 && prediction_bool == 0) {")
  
  # newtxt = write_laplace_estimate(self, private)
  # txt = c(txt, newtxt)
  
  txt = c(txt, "} else if (estimation_method == 3 && prediction_bool == 1) {")
  
  # newtxt = write_laplace_predict(self, private)
  # txt = c(txt, newtxt)
  
  txt = c(txt, "} else {")
  
  txt = c(txt, "}")
  
  # Return nll
  txt = c(txt, "\n //////////// Return //////////////")
  txt = c(txt, "return nll__;")
  txt = c(txt, "}")
  
  # Write cpp function and close file connection
  writeLines(txt,fileconn)
  close(fileconn)
  
  # Return
  return(invisible(self))
}

# txt = c(txt, "\n//////////// UKF METHOD ///////////")
# txt = c(txt, "} else if(estMethod__ == 2){")
# 
# ##################################################
# # UKF OBJECTIVE FUNCTION
# ##################################################
# 
# ##################################################
# # Observation Vectors
# for(i in 1:length(private$obs.names)){
#   txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$obs.names[i]))
# }
# # Input Vectors
# for(i in 1:length(private$input.names)){
#   txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$input.names[i]))
# }
# # Initialize State
# txt = c(txt, "\t DATA_VECTOR(X0__);")
# txt = c(txt, "\t DATA_MATRIX(P0__);")
# # Time-step
# txt = c(txt, "\t DATA_VECTOR(dt__);")
# txt = c(txt, "\t DATA_IVECTOR(N__);")
# # Loss parameters
# txt = c(txt, "\t DATA_VECTOR(tukey_loss_parameters);")
# txt = c(txt, "\t DATA_INTEGER(which_loss__);")
# txt = c(txt, "\t DATA_SCALAR(loss_c_value__);")
# # Maximum a Posterior
# txt = c(txt, "\t DATA_INTEGER(map_bool__);")
# 
# ##################################################
# # Parameters
# for(i in 1:length(private$parameters)){
#   txt = c(txt, sprintf("\t PARAMETER(%s);",private$parameter.names[i]))
# }
# 
# ##################################################
# # Misc
# txt = c(txt, sprintf("\n\t int n__ = %s;",private$n))
# txt = c(txt, sprintf("\t int nn__ = %s;",2*private$n+1))
# txt = c(txt, sprintf("\t int m__ = %s;",private$m))
# txt = c(txt, "\t int s__;")
# 
# ##################################################
# # Storage variables
# txt = c(txt, "\n\t vector<vector<Type>> xPrior(t.size());")
# txt = c(txt, "\t vector<matrix<Type>> pPrior(t.size());")
# txt = c(txt, "\t vector<vector<Type>> xPost(t.size());")
# txt = c(txt, "\t vector<matrix<Type>> pPost(t.size());")
# txt = c(txt, "\t vector<vector<Type>> Innovation(t.size());")
# txt = c(txt, "\t vector<matrix<Type>> InnovationCovariance(t.size());")
# 
# txt = c(txt, "\t xPrior(0) = X0__;")
# txt = c(txt, "\t pPrior(0) = P0__;")
# txt = c(txt, "\t xPost(0) = X0__;")
# txt = c(txt, "\t pPost(0) = P0__;")
# 
# ##################################################
# # Initiaze variables
# txt = c(txt, "\n\t vector<Type> data_vector__(m__);")
# txt = c(txt, "\t vector<Type> na_bool__,y__,e__,X1__;")
# txt = c(txt, "\t matrix<Type> Xsp0__,S0__,G__, F__, S0Inv__, M__, S0PhiM__, Frhs1__(n__,nn__), Frhs0__, Frhs__;")
# txt = c(txt, "\t matrix<Type> E__, H__, Syy__, SyyInv__, Sxy__, K__, P1__;")
# 
# 
# ##################################################
# # Observation equation
# txt = c(txt,"\t matrix<Type> V__(m__,m__);")
# txt = c(txt,"\t V__.setZero();")
# obs.lhs = paste(unlist(lapply(private$obs.eqs.trans,function(x) deparse(x$form[[2]]))),"(i+1)",sep="")
# 
# ##################################################
# # Weights
# txt = c(txt, "\n\t Type lambda__ = 3 - n__, weights__;")
# txt = c(txt, "\t vector<Type> Wm__(nn__);")
# txt = c(txt, "\t matrix<Type> Wmm__(nn__,nn__), Wm_diag__(nn__,nn__), I__(nn__,nn__), W__;")
# txt = c(txt, "\t weights__ = Type(1.0)/(Type(2.0)*(lambda__ + n__));")
# txt = c(txt, "\t Wm__.fill(weights__);")
# txt = c(txt, "\t Wm__(0) = lambda__/(lambda__ + n__);")
# txt = c(txt, "\t for(int i=0; i<nn__ ; i++){")
# txt = c(txt, "\t\t Wmm__.col(i) = Wm__;")
# txt = c(txt, "\t }")
# txt = c(txt, "\t Wm_diag__.setZero();")
# txt = c(txt, "\t Wm_diag__.diagonal() = Wm__;")
# txt = c(txt, "\t I__.setIdentity();")
# txt = c(txt, "\t W__ = (I__ - Wmm__) * Wm_diag__ * (I__ - Wmm__).transpose();\n")
# 
# txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
# txt = c(txt, "\t //////////// MAIN FOR-LOOP ///////////")
# txt = c(txt, "\t ////////////////////////////////////////////////////////")
# 
# txt = c(txt, "\n\t S0__ = P0__.llt().matrixL();")
# txt = c(txt, "\t Xsp0__ = construct_Xsp__(X0__,S0__);")
# 
# ##################################################
# # Time for-loop
# txt = c(txt, "\n\t for(int i=0 ; i<t.size()-1 ; i++){")
# 
# ##################################################
# # Solve Moment ODEs
# txt = c(txt, "\n\t\t //////////// Time-Update ///////////")
# txt = c(txt, "\t\t for(int j=0 ; j<N__(i) ; j++){")
# txt = c(txt, sprintf("\t\t\t F__  = construct_F__(Xsp0__,%s);",paste(fvars2_sigma2,collapse=", ")))
# txt = c(txt, sprintf("\t\t\t G__  = g__(%s);",paste(gvars2_sigma,collapse=", ")) )
# txt = c(txt, "\t\t\t S0Inv__ = S0__.inverse();")
# txt = c(txt, "\t\t\t M__ = S0Inv__ * (Xsp0__ * W__ * F__.transpose() + F__ * W__ * Xsp0__.transpose() + G__*G__.transpose()) * S0Inv__.transpose();")
# txt = c(txt, "\t\t\t S0PhiM__ = S0__ * Phi__(M__);")
# txt = c(txt, "\t\t\t Frhs1__.block(0,1,n__,n__) = S0PhiM__;")
# txt = c(txt, "\t\t\t Frhs1__.block(0,n__+1,n__,n__) = -S0PhiM__;")
# txt = c(txt, "\t\t\t Frhs0__ = (F__*Wm__).replicate(1,nn__);")
# txt = c(txt, "\t\t\t Frhs__ = Frhs0__ + sqrt(3.0) * Frhs1__;")
# txt = c(txt, "\t\t\t Xsp0__ += Frhs__ * dt__(i);")
# txt = c(txt, "\t\t\t S0__ = ((Xsp0__ - Xsp0__.col(0).replicate(1,nn__))/sqrt(Type(3.0))).block(0,1,n__,n__);")
# txt = c(txt, "\t\t\t};")
# txt = c(txt, "\t\t P1__ = S0__ * S0__.transpose();")
# txt = c(txt, "\t\t xPrior(i+1) = Xsp0__.col(0);;")
# txt = c(txt, "\t\t pPrior(i+1) = P1__;")
# 
# txt = c(txt, "\n\t\t //////////// Time-Update ///////////")
# txt = c(txt, sprintf("\t\t data_vector__ << %s;", paste(obs.lhs,collapse=", ")))
# txt = c(txt, "\t\t na_bool__ = is_not_na(data_vector__);")
# txt = c(txt, "\t\t s__ = CppAD::Integer(sum(na_bool__));")
# txt = c(txt, "\t\t if( s__ > 0 ){")
# txt = c(txt, "\t\t\t y__ = remove_nas__(data_vector__,s__,na_bool__);")
# txt = c(txt, "\t\t\t E__ 	= construct_permutation_matrix(s__,m__,na_bool__);")
# txt = c(txt, sprintf("\t\t\t H__ = construct_H__(Xsp0__,%s);", hvars2.withoutStates))
# txt = c(txt, "\t\t\t e__  = y__ - E__ * (H__ * Wm__);")
# txt = c(txt, sprintf("\t\t\t V__ = obsvarFun_usingXsp__(Xsp0__,%s);", obsvars2.withoutStates))
# # txt = c(txt, sprintf("\t\t\t V__.diagonal() << obsvar_diagonal_usingXsp__(Xsp0__, %s);", obsvars2.withoutStates))
# txt = c(txt, "\t\t\t Syy__  = E__ * (H__ * W__ * H__.transpose() + V__) * E__.transpose();")
# txt = c(txt, "\t\t\t Sxy__  = Xsp0__ * W__ * H__.transpose() * E__.transpose();")
# txt = c(txt, "\t\t\t SyyInv__  = Syy__.inverse();")
# txt = c(txt, "\t\t\t K__ = Sxy__ * SyyInv__;")
# txt = c(txt, "\t\t\t X1__ = Xsp0__ * Wm__ + K__ * e__;")
# txt = c(txt, "\t\t\t P1__ = S0__ * S0__.transpose() - K__ * Syy__ * K__.transpose();")
# txt = c(txt, "\t\t\t S0__ = P1__.llt().matrixL();")
# txt = c(txt, "\t\t\t Xsp0__ = construct_Xsp__(X1__,S0__);")
# txt = c(txt, "\t\t\t nll__ += Type(0.5)*atomic::logdet(Syy__) + 0.5*lossfunction__((e__*(SyyInv__*e__)).sum(),tukey_loss_parameters,loss_c_value__,which_loss__) + Type(0.5)*log(2*M_PI)*asDouble(s__);")
# txt = c(txt, "\t\t\t Innovation(i+1) = e__;")
# txt = c(txt, "\t\t\t InnovationCovariance(i+1) = Syy__;")
# txt = c(txt, "\t\t };")
# txt = c(txt, "\t\t xPost(i+1) = Xsp0__.col(0);")
# txt = c(txt, "\t\t pPost(i+1) = P1__;")
# txt = c(txt, "\t };")
# 
# ##################################################
# # Maximum-A-Posterior
# txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
# txt = c(txt, "\t //////////// MAP ESTIMATE USING PRIOR INFO ///////////")
# txt = c(txt, "\t ////////////////////////////////////////////////////////")
# txt = c(txt, "\t if(map_bool__==1){")
# txt = c(txt, "\t\t DATA_VECTOR(map_mean__);")
# txt = c(txt, "\t\t DATA_MATRIX(map_cov__);")
# txt = c(txt, "\t\t DATA_IVECTOR(map_ints__);")
# txt = c(txt, "\t\t DATA_INTEGER(sum_map_ints__);")
# txt = c(txt, sprintf("\t\t vector<Type> parvec__(%s);",length(private$parameters)))
# txt = c(txt, "\t\t vector<Type> map_pars__;")
# txt = c(txt, sprintf("\t\t parvec__ << %s;",paste(private$parameter.names,collapse=", ")))
# txt = c(txt, sprintf("\t\t map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);"))
# txt = c(txt, "\t\t vector<Type> pars_eps__ = map_pars__ - map_mean__;")
# txt = c(txt, "\t\t matrix<Type> map_invcov__ = map_cov__.inverse();")
# txt = c(txt, "\t\t Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
# txt = c(txt, "\t\t nll__ += map_nll__;")
# txt = c(txt, "\t\t REPORT(map_nll__);")
# txt = c(txt, "\t\t REPORT(map_pars__);")
# txt = c(txt, "\t\t REPORT(pars_eps__);")
# txt = c(txt, "\t }")
# 
# ##################################################
# # Report variables and return nll
# txt = c(txt, "\n\t ////////////////////////////////////////////////////////")
# txt = c(txt, "\t //////////// FINAL REPORTING AND RETURN //////////////")
# txt = c(txt, "\t ////////////////////////////////////////////////////////")
# txt = c(txt ,"\t REPORT(xPrior);")
# txt = c(txt, "\t REPORT(pPrior);")
# txt = c(txt,"\t REPORT(xPost);")
# txt = c(txt,"\t REPORT(pPost);")
# txt = c(txt,"\t REPORT(Innovation);")
# txt = c(txt,"\t REPORT(InnovationCovariance);")


# txt = c(txt, "\n//////////// TMB METHOD ///////////")
# txt = c(txt, "} else if (estMethod__ == 3) {")
# 
# ##################################################
# # TMB OBJECTIVE FUNCTION
# ##################################################
# 
# # Observation Vectors
# txt = c(txt, "\n//// observations ////")
# for(i in 1:length(private$obs.names)){
#   txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$obs.names[i]))
# }
# 
# # Input Vectors
# txt = c(txt, "\n//// inputs ////")
# for(i in 1:length(private$input.names)){
#   txt = c(txt, sprintf("\t DATA_VECTOR(%s);",private$input.names[i]))
# }
# 
# # State Vectors
# txt = c(txt, "\n//// state random effect vectors ////")
# for(i in 1:length(private$state.names)){
#   txt = c(txt, sprintf("\t PARAMETER_VECTOR(%s);",private$state.names[i]))
# }
# 
# # Initialize State
# txt = c(txt, "\n//// initial state ////")
# txt = c(txt, "\t DATA_VECTOR(X0__);")
# txt = c(txt, "\t DATA_MATRIX(P0__);")
# 
# # Time-step
# txt = c(txt, "\n//// time-step ////")
# txt = c(txt, "\t DATA_VECTOR(dt__);")
# txt = c(txt, "\t DATA_IVECTOR(N__);")
# txt = c(txt, "\t DATA_IVECTOR(Nc__);")
# 
# # Iobs
# txt = c(txt, "\n//// iobs vectors ////")
# for (i in 1:private$m) {
#   nam = paste("iobs_",private$obs.names[i],sep="")
#   txt = c(txt, sprintf("\t DATA_IVECTOR(%s);",nam))
# }
# 
# # Loss parameters
# # txt = c(txt, "\t DATA_VECTOR(tukey_loss_parameters);")
# # txt = c(txt, "\t DATA_INTEGER(which_loss__);")
# # txt = c(txt, "\t DATA_SCALAR(loss_c_value__);")
# 
# # Maximum a Posterior
# txt = c(txt, "\n//// map estimation ////")
# txt = c(txt, "\t DATA_INTEGER(map_bool__);")
# 
# # Parameters
# txt = c(txt, "\n//// parameters ////")
# for(i in 1:length(private$parameters)){
#   txt = c(txt, sprintf("\t PARAMETER(%s);",private$parameter.names[i]))
# }
# 
# # system size
# txt = c(txt, "\n//// system size ////")
# txt = c( txt , "\t DATA_INTEGER(n__);")
# txt = c( txt , "\t DATA_INTEGER(m__);")
# txt = c( txt , "\t DATA_INTEGER(ng__);")
# 
# # Likelihood
# # txt = c(txt,sprintf("\n\t int n__ = %s;",private$n))
# # txt = c(txt,sprintf("\t int ng__ = %s;",private$ng))
# # txt = c(txt,sprintf("\t int m__ = %s;",private$m))
# # txt = c(txt,"\t Type nll__ = 0;")
# 
# # Initiaze variables
# txt = c(txt, "\n\t //////////// initialize variables ///////////")
# txt = c(txt, "\t vector<Type> F__(n__);")
# txt = c(txt, "\t vector<Type> Xi__(n__);")
# txt = c(txt, "\t vector<Type> Xip1__(n__);")
# txt = c(txt, "\t vector<Type> Z__(n__);")
# txt = c(txt, "\t matrix<Type> G__(n__,ng__);")
# txt = c(txt, "\t matrix<Type> V__(n__,n__);")
# txt = c(txt, "\t matrix<Type> I__(n__,n__);")
# txt = c(txt, "\t I__.setIdentity();")
# txt = c(txt, "\t I__ *= 1e-8;")
# 
# # Forward simulation and likelihood contribution from states
# txt = c(txt, "\n\t //////////// MAIN LOOP OVER TIME POINTS ///////////")
# 
# txt = c(txt, "\t for(int i=0 ; i<t.size()-1 ; i++){")
# #
# txt = c(txt, "\t\t for(int j=0 ; j<N__(i) ; j++){")
# #
# txt = c(txt, sprintf("\t\t\t F__ = f__(%s);",paste(fvars2.tmb,collapse=", ")))
# txt = c(txt, sprintf("\t\t\t G__ = g__(%s);",paste(gvars2.tmb,collapse=", ")))
# txt = c(txt, sprintf("\t\t\t Xi__ << %s;",paste(private$state.names,"(Nc__(i)+j)",collapse=", ",sep="")))
# txt = c(txt, sprintf("\t\t\t Xip1__ << %s;",paste(private$state.names,"(Nc__(i)+j+1)",collapse=", ",sep="")))
# txt = c(txt, sprintf("\t\t\t Z__ = Xip1__ - Xi__ - F__ * dt__(i);"))
# txt = c(txt, sprintf("\t\t\t V__ = (G__ * G__.transpose() + I__) * dt__(i);"))
# txt = c(txt, "\t\t\t nll__ += MVNORM(V__)(Z__);")
# #
# txt = c(txt, "\t\t }")
# #
# txt = c(txt, "\t }")
# 
# # Data Update
# txt = c(txt, "\n\t\t //////////// DATA-UPDATE ///////////")
# #
# txt = c(txt, "\t matrix<Type> varDiag__(m__,t.size());")
# txt = c(txt, "\t matrix<Type> varFun__(m__,t.size());")
# txt = c(txt, "\t for(int i=0 ; i<t.size() ; i++){")
# txt = c(txt, sprintf("\t\t varDiag__.col(i) = obsvarFun_tmb__(%s);",paste(obsvars2.tmb,collapse=", ")))
# txt = c(txt, sprintf("\t\t varFun__.col(i) = h__(%s);",paste(hvars2.tmb,collapse=", ")))
# txt = c(txt, "\t }")
# #
# for(i in 1:private$m){
#   nam = paste("iobs_",private$obs.names[i],sep="")
#   txt = c(txt, sprintf("\t for(int i=0 ; i<%s.size() ; i++){",nam))
#   txt = c(txt, sprintf("\t\t int j = %s(i);",nam))
#   txt = c(txt, sprintf("\t\t nll__ -= dnorm(%s,varFun__.col(j)(%s),sqrt(varDiag__.col(j)(%s)),true);",
#                        paste(private$obs.names[i],"(j)",sep=""),
#                        i-1,
#                        i-1))
#   txt = c(txt, "\t }")
# }
# 
# # Maximum-A-Posterior
# txt = c(txt, "\n\t\t //////////// MAP CONTRIBUTION ///////////")
# #
# txt = c(txt, "\t if(map_bool__==1){")
# txt = c(txt, "\t\t DATA_VECTOR(map_mean__);")
# txt = c(txt, "\t\t DATA_MATRIX(map_cov__);")
# txt = c(txt, "\t\t DATA_IVECTOR(map_ints__);")
# txt = c(txt, "\t\t DATA_INTEGER(sum_map_ints__);")
# txt = c(txt, sprintf("\t\t vector<Type> parvec__(%s);",length(private$parameters)))
# txt = c(txt, "\t\t vector<Type> map_pars__;")
# txt = c(txt, sprintf("\t\t parvec__ << %s;",paste(private$parameter.names,collapse=", ")))
# txt = c(txt, sprintf("\t\t map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);"))
# txt = c(txt, "\t\t vector<Type> pars_eps__ = map_pars__ - map_mean__;")
# txt = c(txt, "\t\t matrix<Type> map_invcov__ = map_cov__.inverse();")
# txt = c(txt, "\t\t Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();")
# txt = c(txt, "\t\t nll__ += map_nll__;")
# txt = c(txt, "\t\t REPORT(map_nll__);")
# txt = c(txt, "\t\t REPORT(map_pars__);")
# txt = c(txt, "\t\t REPORT(pars_eps__);")
# txt = c(txt, "\t }")
# txt = c(txt, "}")

#   txt = c(txt, "\n //////////// Estimate Methods Ends //////////////")
#   txt = c(txt, "} else {")
#   txt = c(txt, "Type nll = 0;")
#   txt = c(txt, "}")
#   
#   # Return nll
#   txt = c(txt, "\n //////////// Return //////////////")
#   txt = c(txt, "return nll__;")
#   txt = c(txt, "}")
#   
#   # Write cpp function and close file connection
#   writeLines(txt,fileconn)
#   close(fileconn)
#   
#   # Return
#   return(invisible(self))
# }

