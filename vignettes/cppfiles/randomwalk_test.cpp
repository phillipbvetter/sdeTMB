#include <TMB.hpp>
using namespace density;

//////////// helper fun: find NA locations in vector ///////////
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
}

//////////// helper fun: extract non-NAs from vector ///////////
template<class Type>
vector<Type> remove_nas__(vector<Type> data_vector, int number_of_datas, vector<Type> na_bool){
  int ii = 0;
	vector<Type> y_reduced(number_of_datas);
	for(int i=0; i < data_vector.size(); i++){
		if(na_bool(i) == Type(1.0)){
			y_reduced(ii) = data_vector(i);
			ii++;
		}
	}
  return y_reduced;
}

//////////// helper fun: construct permutation matrix ///////////
template <class Type>
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
}

//////////// loss function ///////////
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
  return(loss);
}

//////////// MAP estimation helper ///////////
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
	return(ans);
}

//////////// UKF sigma points ///////////
template<class Type>
matrix<Type> construct_Xsp__(vector<Type> x0, matrix<Type> s0){
  int n = 1;
  int nn = 3;
  matrix<Type> Xsp(n,nn);
  vector<Type> Si;
  Xsp.col(0) = x0;
  for(int i=1; i<n+1; i++){
    Si = s0.col(i-1);
    Xsp.col(i) = x0 + sqrt(1+n) * Si;
    Xsp.col(i+n) = x0 - sqrt(1+n) * Si;
  }
  return Xsp;
}

//////////// UKF function ///////////
template<class Type>
matrix<Type> Phi__(matrix<Type> M){
  matrix<Type> K(M.col(0).size(),M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/Type(2.0);
  return K;
}

//////////// drift function ///////////
template<class Type>
vector<Type> f__(Type a){
	vector<Type> ans(1);
	ans(0) = 0;
	return ans;
}

//////////// UKF sigma points drift function ///////////
template<class Type>
matrix<Type> construct_F__(matrix<Type> Xsp, Type a){
  int n = 1;
  int nn = 3;
  matrix<Type> F(n,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    F.col(i) = f__(Type(0.0));
  }
  return F;
}

//////////// jacobian of drift function ///////////
template<class Type>
matrix<Type> dfdx__(Type void_filler){
	matrix<Type> ans(1,1);
	ans(0,0) = 0;
	return ans;
}

//////////// diffusion function ///////////
template<class Type>
matrix<Type> g__(Type logsigma_x){
	matrix<Type> ans(1,1);
	ans(0,0) = exp(logsigma_x);
	return ans;
}

//////////// observation function ///////////
template<class Type>
vector<Type> h__(Type x){
	vector<Type> ans(1);
	ans(0) = x;
	return ans;
}

//////////// jacobian of observation function ///////////
template<class Type>
matrix<Type> dhdx__(Type state_void_filler){
	matrix<Type> ans(1,1);
	ans(0,0) = 1;
	return ans;
}

//////////// UKF sigma points observation function ///////////
template<class Type>
matrix<Type> construct_H__(matrix<Type> Xsp, Type void_filler){
  int nn = 3;
  int m = 1;
  matrix<Type> H(m,nn);
  vector<Type> x0__;
  for(int i=0;i<nn;i++){
  x0__ = Xsp.col(i);
  H.col(i) = h__(x0__(0));
  }
return H;
}

//////////// observation variance matrix function ///////////
template<class Type>
matrix<Type> obsvarFun__(Type logsigma_y){
	matrix<Type> V(1,1);
	V(0,0) = exp(pow(logsigma_y, 2));
	return V;
}

//////////// observation variance vector function ///////////
template<class Type>
vector<Type> obsvarFun_tmb__(Type logsigma_y){
	vector<Type> ans(1);
	ans(0) = exp(pow(logsigma_y, 2));
	return ans;
}

//////////// UKF sigma points observation variance matrix function ///////////
template<class Type>
matrix<Type> obsvarFun_usingXsp__(matrix<Type> Xsp, Type logsigma_y){
  matrix<Type> V;
  vector<Type> x0__ = Xsp.col(0);
  V = obsvarFun__(logsigma_y);
  return V;
  }

//////////// objective function ///////////
template<class Type>
Type objective_function<Type>::operator() ()
{
	 DATA_INTEGER(estMethod__);
	 DATA_INTEGER(pred__);
	 Type nll__ = 0;

//////////// EKF METHOD ///////////
	 if(estMethod__ == 1){

	 //////////// Estimation //////////////

	 //////////// Estimation //////////////
	 if(pred__ == 0){

//// observations ////
	 DATA_VECTOR(y);

//// inputs ////
	 DATA_VECTOR(t);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);

//// loss parameters ////
	 DATA_VECTOR(tukey_pars__);
	 DATA_INTEGER(which_loss__);
	 DATA_SCALAR(loss_c_value__);

//// map estimation ////
	 DATA_INTEGER(map_bool__);

//// parameters ////
	 PARAMETER(logsigma_x);
	 PARAMETER(logsigma_y);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);

//////////// storage variables ///////////
	 vector<vector<Type>> xPrior(t.size());
	 vector<matrix<Type>> pPrior(t.size());
	 vector<vector<Type>> xPost(t.size());
	 vector<matrix<Type>> pPost(t.size());
	 vector<vector<Type>> Innovation(t.size());
	 vector<matrix<Type>> InnovationCovariance(t.size());

//////////// set initial value ///////////
	 vector<Type> x0__ = X0__;
	 matrix<Type> p0__ = P0__;
	 xPrior(0) = X0__;
	 xPost(0) = X0__;
	 pPrior(0) = P0__;
	 pPost(0) = P0__;

	 //////////// initialize variables ///////////
	 int s__;
	 Type half_log2PI = Type(0.5)*log(2*M_PI);
	 vector<Type> data_vector__(m__),na_bool__,e__,y__,F__,H__;
	 matrix<Type> C__,R__,K__,E__,V__,Ri__,A__,G__;

	 //////////// identity matrix ///////////
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<t.size()-1 ; i++){

		 //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
		 for(int j=0 ; j<N__(i) ; j++){
			 F__  = f__(Type(0.0));
			 A__  = dfdx__(Type(0.0));
			 G__  = g__(logsigma_x);
			 x0__ = x0__ + F__ * dt__(i);
			 p0__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt__(i);
		 }
		 xPrior(i+1) = x0__;
		 pPrior(i+1) = p0__;

		 //////////// DATA-UPDATE ///////////
		 data_vector__ << y(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas__(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(logsigma_y) * E__.transpose();
			 R__  = C__ * p0__ * C__.transpose() + V__;
			 Ri__ = R__.inverse();
			 K__ 	= p0__ * C__.transpose() * Ri__;
			 x0__ = x0__ + K__*e__;
			 p0__ = (I__ - K__ * C__) * p0__ * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();
			 nll__ += Type(0.5)*atomic::logdet(R__) + Type(0.5)*lossfunction__((e__*(Ri__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + half_log2PI * asDouble(s__);
			 Innovation(i+1) = e__;
			 InnovationCovariance(i+1) = R__;
		 }
		 xPost(i+1) = x0__;
		 pPost(i+1) = p0__;
	 }

		 //////////// MAP CONTRIBUTION ///////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(2);
		 vector<Type> map_pars__;
		 parvec__ << logsigma_x, logsigma_y;
		 map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);
		 vector<Type> pars_eps__ = map_pars__ - map_mean__;
		 matrix<Type> map_invcov__ = map_cov__.inverse();
		 Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
		 nll__ += map_nll__;
		 REPORT(map_nll__);
		 REPORT(map_pars__);
		 REPORT(pars_eps__);
	 }

	 //////////// Return/Report //////////////
	 REPORT(Innovation);
	 REPORT(InnovationCovariance);
	 REPORT(xPrior);
	 REPORT(xPost);
	 REPORT(pPrior);
	 REPORT(pPost);

	 //////////// Prediction //////////////

	 //////////// Prediction //////////////
	 } else if(pred__ == 1){

//// observations ////
	 DATA_VECTOR(y);

//// inputs ////
	 DATA_VECTOR(t);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);

//// parameters ////
	 PARAMETER(logsigma_x);
	 PARAMETER(logsigma_y);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);
	 DATA_INTEGER(last_pred_index);
	 DATA_INTEGER(k_step_ahead);
	 vector<matrix<Type>> xk__(last_pred_index);
	 vector<matrix<Type>> pk__(last_pred_index);
	 xk__.fill(matrix<Type>(k_step_ahead+1,n__));
	 pk__.fill(matrix<Type>(k_step_ahead+1,n__*n__));
	 matrix<Type> xk_temp__(k_step_ahead+1,n__);
	 array<Type> pk_temp__(n__,n__,k_step_ahead+1);

//////////// set initial value ///////////
	 vector<Type> x0__ = X0__;
	 matrix<Type> p0__ = P0__;

	 //////////// initialize variables ///////////
	 int s__;
	 vector<Type> data_vector__(m__),na_bool__,e__,y__,F__,H__;
	 matrix<Type> C__,R__,K__,E__,V__,Ri__,A__,G__;

	 //////////// identity matrix ///////////
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<last_pred_index ; i++){

		 xk_temp__.row(0) = x0__;
		 pk_temp__.col(0) = p0__;

	 //////////// K-STEP-AHEAD LOOP ///////////
	 for(int k=0 ; k < k_step_ahead ; k++){

		 //////////// TIME-UPDATE: SOLVE MOMENT ODES ///////////
		 for(int j=0 ; j<N__(i+k) ; j++){
			 F__  = f__();
			 A__  = dfdx__(Type(0.0));
			 G__  = g__(logsigma_x);
			 x0__ = x0__ + F__ * dt__(i+k);
			 p0__ = p0__ + ( A__*p0__ + p0__*A__.transpose() + G__*G__.transpose() ) * dt__(i+k);
		 }

			 //////////// save k-step-ahead prediction ///////////
			 xk_temp__.row(k+1) = x0__;
			 pk_temp__.col(k+1) = p0__;
		 }

		 //////////// save all 0 to k step-ahead predictions ///////////
		 xk__(i) = xk_temp__;
		 for(int kk=0 ; kk < k_step_ahead+1 ; kk++){
			 pk__(i).row(kk) = vector<Type>(pk_temp__.col(kk).transpose());
		 }

		 //////////// rewrite x0 and p0 to one-step predictions ///////////
		 x0__ = xk_temp__.row(1);
		 p0__ = pk_temp__.col(1).matrix();

		 //////////// DATA-UPDATE ///////////
		 data_vector__ << y(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas__(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(logsigma_y) * E__.transpose();
			 R__  = C__ * p0__ * C__.transpose() + V__;
			 Ri__ = R__.inverse();
			 K__ 	= p0__ * C__.transpose() * Ri__;
			 x0__ = x0__ + K__*e__;
			 p0__ = (I__ - K__ * C__) * p0__ * (I__ - K__ * C__).transpose() + K__* V__ * K__.transpose();
		 }
	 }

	 //////////// Return/Report //////////////
	 REPORT(xk__);
	 REPORT(pk__);
	 }

//////////// UKF METHOD ///////////
} else if(estMethod__ == 2){
	 DATA_VECTOR(y);
	 DATA_VECTOR(t);
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);
	 DATA_VECTOR(tukey_pars__);
	 DATA_INTEGER(which_loss__);
	 DATA_SCALAR(loss_c_value__);
	 DATA_INTEGER(map_bool__);
	 PARAMETER(logsigma_x);
	 PARAMETER(logsigma_y);

	 int n__ = 1;
	 int nn__ = 3;
	 int m__ = 1;
	 int s__;

	 vector<vector<Type>> xPrior(t.size());
	 vector<matrix<Type>> pPrior(t.size());
	 vector<vector<Type>> xPost(t.size());
	 vector<matrix<Type>> pPost(t.size());
	 vector<vector<Type>> Innovation(t.size());
	 vector<matrix<Type>> InnovationCovariance(t.size());
	 xPrior(0) = X0__;
	 pPrior(0) = P0__;
	 xPost(0) = X0__;
	 pPost(0) = P0__;

	 vector<Type> data_vector__(m__);
	 vector<Type> na_bool__,y__,e__,X1__;
	 matrix<Type> Xsp0__,S0__,G__, F__, S0Inv__, M__, S0PhiM__, Frhs1__(n__,nn__), Frhs0__, Frhs__;
	 matrix<Type> E__, H__, Syy__, SyyInv__, Sxy__, K__, P1__;
	 matrix<Type> V__(m__,m__);
	 V__.setZero();

	 Type lambda__ = 3 - n__, weights__;
	 vector<Type> Wm__(nn__);
	 matrix<Type> Wmm__(nn__,nn__), Wm_diag__(nn__,nn__), I__(nn__,nn__), W__;
	 weights__ = Type(1.0)/(Type(2.0)*(lambda__ + n__));
	 Wm__.fill(weights__);
	 Wm__(0) = lambda__/(lambda__ + n__);
	 for(int i=0; i<nn__ ; i++){
		 Wmm__.col(i) = Wm__;
	 }
	 Wm_diag__.setZero();
	 Wm_diag__.diagonal() = Wm__;
	 I__.setIdentity();
	 W__ = (I__ - Wmm__) * Wm_diag__ * (I__ - Wmm__).transpose();


	 ////////////////////////////////////////////////////////
	 //////////// MAIN FOR-LOOP ///////////
	 ////////////////////////////////////////////////////////

	 S0__ = P0__.llt().matrixL();
	 Xsp0__ = construct_Xsp__(X0__,S0__);

	 for(int i=0 ; i<t.size()-1 ; i++){

		 //////////// Time-Update ///////////
		 for(int j=0 ; j<N__(i) ; j++){
			 F__  = construct_F__(Xsp0__,Type(0.0));
			 G__  = g__(logsigma_x);
			 S0Inv__ = S0__.inverse();
			 M__ = S0Inv__ * (Xsp0__ * W__ * F__.transpose() + F__ * W__ * Xsp0__.transpose() + G__*G__.transpose()) * S0Inv__.transpose();
			 S0PhiM__ = S0__ * Phi__(M__);
			 Frhs1__.block(0,1,n__,n__) = S0PhiM__;
			 Frhs1__.block(0,n__+1,n__,n__) = -S0PhiM__;
			 Frhs0__ = (F__*Wm__).replicate(1,nn__);
			 Frhs__ = Frhs0__ + sqrt(3.0) * Frhs1__;
			 Xsp0__ += Frhs__ * dt__(i);
			 S0__ = ((Xsp0__ - Xsp0__.col(0).replicate(1,nn__))/sqrt(Type(3.0))).block(0,1,n__,n__);
			};
		 P1__ = S0__ * S0__.transpose();
		 xPrior(i+1) = Xsp0__.col(0);;
		 pPrior(i+1) = P1__;

		 //////////// Time-Update ///////////
		 data_vector__ << y(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__ = remove_nas__(data_vector__,s__,na_bool__);
			 E__ 	= construct_permutation_matrix(s__,m__,na_bool__);
			 H__ = construct_H__(Xsp0__,Type(0.0));
			 e__  = y__ - E__ * (H__ * Wm__);
			 V__ = obsvarFun_usingXsp__(Xsp0__,logsigma_y);
			 Syy__  = E__ * (H__ * W__ * H__.transpose() + V__) * E__.transpose();
			 Sxy__  = Xsp0__ * W__ * H__.transpose() * E__.transpose();
			 SyyInv__  = Syy__.inverse();
			 K__ = Sxy__ * SyyInv__;
			 X1__ = Xsp0__ * Wm__ + K__ * e__;
			 P1__ = S0__ * S0__.transpose() - K__ * Syy__ * K__.transpose();
			 S0__ = P1__.llt().matrixL();
			 Xsp0__ = construct_Xsp__(X1__,S0__);
			 nll__ += Type(0.5)*atomic::logdet(Syy__) + 0.5*lossfunction__((e__*(SyyInv__*e__)).sum(),tukey_pars__,loss_c_value__,which_loss__) + Type(0.5)*log(2*M_PI)*asDouble(s__);
			 Innovation(i+1) = e__;
			 InnovationCovariance(i+1) = Syy__;
		 };
		 xPost(i+1) = Xsp0__.col(0);
		 pPost(i+1) = P1__;
	 };

	 ////////////////////////////////////////////////////////
	 //////////// MAP ESTIMATE USING PRIOR INFO ///////////
	 ////////////////////////////////////////////////////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(2);
		 vector<Type> map_pars__;
		 parvec__ << logsigma_x, logsigma_y;
		 map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);
		 vector<Type> pars_eps__ = map_pars__ - map_mean__;
		 matrix<Type> map_invcov__ = map_cov__.inverse();
		 Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
		 nll__ += map_nll__;
		 REPORT(map_nll__);
		 REPORT(map_pars__);
		 REPORT(pars_eps__);
	 }

	 ////////////////////////////////////////////////////////
	 //////////// FINAL REPORTING AND RETURN //////////////
	 ////////////////////////////////////////////////////////
	 REPORT(xPrior);
	 REPORT(pPrior);
	 REPORT(xPost);
	 REPORT(pPost);
	 REPORT(Innovation);
	 REPORT(InnovationCovariance);

//////////// TMB METHOD ///////////
} else if(estMethod__ == 3) {

//// observations ////
	 DATA_VECTOR(y);

//// inputs ////
	 DATA_VECTOR(t);

//// state random effect vectors ////
	 PARAMETER_VECTOR(x);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);

//// time-step ////
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);
	 DATA_IVECTOR(Nc__);

//// iobs vectors ////
	 DATA_IVECTOR(iobs_y);

//// map estimation ////
	 DATA_INTEGER(map_bool__);

//// parameters ////
	 PARAMETER(logsigma_x);
	 PARAMETER(logsigma_y);

//// constants ////

//// system size ////
	 DATA_INTEGER(n__);
	 DATA_INTEGER(m__);
	 DATA_INTEGER(ng__);

	 //////////// initialize variables ///////////
	 vector<Type> F__(n__);
	 vector<Type> Xi__(n__);
	 vector<Type> Xip1__(n__);
	 vector<Type> Z__(n__);
	 matrix<Type> G__(n__,ng__);
	 matrix<Type> V__(n__,n__);
	 matrix<Type> I__(n__,n__);
	 I__.setIdentity();
	 I__ *= 1e-8;

	 //////////// MAIN LOOP OVER TIME POINTS ///////////
	 for(int i=0 ; i<t.size()-1 ; i++){
		 for(int j=0 ; j<N__(i) ; j++){
			 F__ = f__(Type(0.0));
			 G__ = g__(logsigma_x);
			 Xi__ << x(Nc__(i)+j);
			 Xip1__ << x(Nc__(i)+j+1);
			 Z__ = Xip1__ - Xi__ - F__ * dt__(i);
			 V__ = (G__ * G__.transpose() + I__) * dt__(i);
			 nll__ += MVNORM(V__)(Z__);
		 }
	 }

		 //////////// DATA-UPDATE ///////////
	 matrix<Type> varDiag__(m__,t.size());
	 for(int i=0 ; i<t.size() ; i++){
		 varDiag__.col(i) = obsvarFun_tmb__(logsigma_y);
	 }
	 for(int i=0 ; i<iobs_y.size() ; i++){
		 int j = iobs_y(i);
		 nll__ -= dnorm(y(j),x(Nc__(j)),sqrt(varDiag__.col(j)(0)),true);
	 }

		 //////////// MAP CONTRIBUTION ///////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(2);
		 vector<Type> map_pars__;
		 parvec__ << logsigma_x, logsigma_y;
		 map_pars__ = get_free_pars__(map_ints__,sum_map_ints__,parvec__);
		 vector<Type> pars_eps__ = map_pars__ - map_mean__;
		 matrix<Type> map_invcov__ = map_cov__.inverse();
		 Type map_nll__ = Type(0.5) * atomic::logdet(map_cov__) + Type(0.5) * (pars_eps__ * (map_invcov__ * pars_eps__)).sum();
		 nll__ += map_nll__;
		 REPORT(map_nll__);
		 REPORT(map_pars__);
		 REPORT(pars_eps__);
	 }
}
return nll__;
}
