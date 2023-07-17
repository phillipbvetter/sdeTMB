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
vector<Type> f__(Type mu, Type theta, Type X){
	vector<Type> ans(1);
	ans(0) = theta * (mu - X);
	return ans;
}

//////////// UKF sigma points drift function ///////////
template<class Type>
matrix<Type> construct_F__(matrix<Type> Xsp, Type mu, Type theta){
  int n = 1;
  int nn = 3;
  matrix<Type> F(n,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    F.col(i) = f__(mu, theta, x0(0));
  }
  return F;
}

//////////// jacobian of drift function ///////////
template<class Type>
matrix<Type> dfdx__(Type theta){
	matrix<Type> ans(1,1);
	ans(0,0) = -theta;
	return ans;
}

//////////// diffusion function ///////////
template<class Type>
matrix<Type> g__(Type sigma_x){
	matrix<Type> ans(1,1);
	ans(0,0) = sigma_x;
	return ans;
}

//////////// observation function ///////////
template<class Type>
vector<Type> h__(Type X){
	vector<Type> ans(1);
	ans(0) = X;
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
struct ode_integration {
	vector<Type> X_next;
	matrix<Type> P_next;
	ode_integration(vector<Type> x0__, matrix<Type> p0__, Type t, Type dt__, int algo__, Type mu, Type sigma_x, Type theta){
		if(algo__==1){
			/*Forward Euler*/
			X_next = x0__ + f__(mu, theta, x0__(0)) * dt__;
			P_next = p0__ + (dfdx__(theta)*p0__ + p0__*dfdx__(theta).transpose() + g__(sigma_x)*g__(sigma_x).transpose()) * dt__;
		} else if (algo__==2){
			/*4th Order Runge-Kutta 4th*/
			vector<Type> X0__ = x0__;
			matrix<Type> P0__ = p0__;
			/**/
			vector<Type> k1__,k2__,k3__,k4__;
			matrix<Type> a1__,a2__,a3__,a4__;
			/*SOLVE ODE*/
			/*step 1*/
			k1__ = dt__ * f__(mu, theta, x0__(0));
			a1__ = dt__ * (dfdx__(theta)*p0__ + p0__*dfdx__(theta).transpose() + g__(sigma_x)*g__(sigma_x).transpose());
			/*step 2*/
			t = t + 0.5 * dt__;
			x0__ = X0__ + 0.5 * k1__;
			p0__ = P0__ + 0.5 * a1__;
      k2__ = dt__ * f__(mu, theta, x0__(0));
      a2__ = dt__ * (dfdx__(theta)*p0__ + p0__*dfdx__(theta).transpose() + g__(sigma_x)*g__(sigma_x).transpose());
			/*step 3*/
			x0__ = X0__ + 0.5 * k2__;
			p0__ = P0__ + 0.5 * a2__;
      k3__ = dt__ * f__(mu, theta, x0__(0));
      a3__ = dt__ * (dfdx__(theta)*p0__ + p0__*dfdx__(theta).transpose() + g__(sigma_x)*g__(sigma_x).transpose());
			/*step 4*/
			t = t + 0.5 * dt__;
			x0__ = X0__ + k3__;
			p0__ = P0__ + a3__;
      k4__ = dt__ * f__(mu, theta, x0__(0));
      a4__ = dt__ * (dfdx__(theta)*p0__ + p0__*dfdx__(theta).transpose() + g__(sigma_x)*g__(sigma_x).transpose());
			/*ODE UPDATE*/
			X_next = X0__ + (k1__ + 2.0*k2__ + 2.0*k3__ + k4__)/6.0;
			P_next = P0__ + (a1__ + 2.0*a2__ + 2.0*a3__ + a4__)/6.0;
		} else {
			/*nothing*/
		}
	}
};
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
matrix<Type> obsvarFun__(Type sigma_y, Type U){
	matrix<Type> V(1,1);
	V(0,0) = pow(sigma_y, 2) * U;
	return V;
}

//////////// observation variance vector function ///////////
template<class Type>
vector<Type> obsvarFun_tmb__(Type sigma_y, Type U){
	vector<Type> ans(1);
	ans(0) = pow(sigma_y, 2) * U;
	return ans;
}

//////////// UKF sigma points observation variance matrix function ///////////
template<class Type>
matrix<Type> obsvarFun_usingXsp__(matrix<Type> Xsp, Type sigma_y, Type U){
  matrix<Type> V;
  vector<Type> x0__ = Xsp.col(0);
  V = obsvarFun__(sigma_y, U);
  return V;
  }

//////////// objective function ///////////
template<class Type>
Type objective_function<Type>::operator() ()
{
	 DATA_INTEGER(estMethod__);
	 DATA_INTEGER(pred__);
	 DATA_INTEGER(algo__);
	 Type nll__ = 0;

//////////// EKF METHOD ///////////
	 if(estMethod__ == 1){

	 //////////// Estimation //////////////

	 //////////// Estimation //////////////
	 if(pred__ == 0){

//// observations ////
	 DATA_VECTOR(Y);

//// inputs ////
	 DATA_VECTOR(t);
	 DATA_VECTOR(U);

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
	 PARAMETER(theta);
	 PARAMETER(mu);
	 PARAMETER(sigma_x);
	 PARAMETER(sigma_y);

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
			 ode_integration<Type> odelist = {x0__, p0__, t(i)+j*dt__(i), dt__(i), algo__, mu, sigma_x, theta};
			 x0__ = odelist.X_next;
			 p0__ = odelist.P_next;
		 }
		 xPrior(i+1) = x0__;
		 pPrior(i+1) = p0__;

		 //////////// DATA-UPDATE ///////////
		 data_vector__ << Y(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas__(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(sigma_y ,U(i)) * E__.transpose();
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
		 vector<Type> parvec__(4);
		 vector<Type> map_pars__;
		 parvec__ << theta, mu, sigma_x, sigma_y;
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
	 DATA_VECTOR(Y);

//// inputs ////
	 DATA_VECTOR(t);
	 DATA_VECTOR(U);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);

//// parameters ////
	 PARAMETER(theta);
	 PARAMETER(mu);
	 PARAMETER(sigma_x);
	 PARAMETER(sigma_y);

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
			 ode_integration<Type> odelist = {x0__, p0__, t(i+k)+j*dt__(i+k), dt__(i+k), algo__, mu, sigma_x, theta};
			 x0__ = odelist.X_next;
			 p0__ = odelist.P_next;
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
		 data_vector__ << Y(i+1);
		 na_bool__ = is_not_na(data_vector__);
		 s__ = CppAD::Integer(sum(na_bool__));
		 if( s__ > 0 ){
			 y__  = remove_nas__(data_vector__, s__, na_bool__);
			 E__  = construct_permutation_matrix(s__, m__, na_bool__);
			 H__  = h__(x0__(0));
			 C__  = E__ * dhdx__(Type(0.0));
			 e__  = y__ - E__ * H__;
			 V__  = E__ * obsvarFun__(sigma_y ,U(i)) * E__.transpose();
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

//////////// TMB METHOD ///////////
} else if(estMethod__ == 3) {

//// observations ////
	 DATA_VECTOR(Y);

//// inputs ////
	 DATA_VECTOR(t);
	 DATA_VECTOR(U);

//// state random effect vectors ////
	 PARAMETER_VECTOR(X);

//// initial state ////
	 DATA_VECTOR(X0__);
	 DATA_MATRIX(P0__);

//// time-step ////
	 DATA_VECTOR(dt__);
	 DATA_IVECTOR(N__);
	 DATA_IVECTOR(Nc__);

//// iobs vectors ////
	 DATA_IVECTOR(iobs_Y);

//// map estimation ////
	 DATA_INTEGER(map_bool__);

//// parameters ////
	 PARAMETER(theta);
	 PARAMETER(mu);
	 PARAMETER(sigma_x);
	 PARAMETER(sigma_y);

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
			 F__ = f__(mu, theta, X(Nc__(i)+j));
			 G__ = g__(sigma_x);
			 Xi__ << X(Nc__(i)+j);
			 Xip1__ << X(Nc__(i)+j+1);
			 Z__ = Xip1__ - Xi__ - F__ * dt__(i);
			 V__ = (G__ * G__.transpose() + I__) * dt__(i);
			 nll__ += MVNORM(V__)(Z__);
		 }
	 }

		 //////////// DATA-UPDATE ///////////
	 matrix<Type> varDiag__(m__,t.size());
	 matrix<Type> varFun__(m__,t.size());
	 for(int i=0 ; i<t.size() ; i++){
		 varDiag__.col(i) = obsvarFun_tmb__(sigma_y, U(i));
		 varFun__.col(i) = h__(X(Nc__(i)));
	 }
	 for(int i=0 ; i<iobs_Y.size() ; i++){
		 int j = iobs_Y(i);
		 nll__ -= dnorm(Y(j),varFun__.col(j)(0),sqrt(varDiag__.col(j)(0)),true);
	 }

		 //////////// MAP CONTRIBUTION ///////////
	 if(map_bool__==1){
		 DATA_VECTOR(map_mean__);
		 DATA_MATRIX(map_cov__);
		 DATA_IVECTOR(map_ints__);
		 DATA_INTEGER(sum_map_ints__);
		 vector<Type> parvec__(4);
		 vector<Type> map_pars__;
		 parvec__ << theta, mu, sigma_x, sigma_y;
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
