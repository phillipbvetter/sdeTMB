#include <TMB.hpp>
#include <cmath>
using namespace density;
template <class Type>
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
}
template<class Type>
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
}
template <class Type>
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
}
template<class Type>
matrix<Type> construct_Xsp(vector<Type> x0, matrix<Type> s0){
  int n = 2;
  int nn = 5;
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
template<class Type>
matrix<Type> Phi(matrix<Type> M){
  matrix<Type> K(M.col(0).size(),M.row(0).size());
  K.setZero();
  K.template triangularView<Eigen::Lower>() = M.template triangularView<Eigen::Lower>();
  K.diagonal() = K.diagonal()/Type(2.0);
  return K;
}
template<class Type>
vector<Type> f(Type loglambda, Type x1, Type x2){
	vector<Type> ans(2);
	ans(0) = x2;
	ans(1) = x2 * (1 - pow(x1, 2)) * exp(loglambda) - x1;
	return ans;
}
template<class Type>
matrix<Type> construct_F(matrix<Type> Xsp, Type loglambda){
  int n = 2;
  int nn = 5;
  matrix<Type> F(n,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    F.col(i) = f(loglambda, x0(0), x0(1));
  }
  return F;
}
template<class Type>
matrix<Type> g(Type logsigma_x1, Type logsigma_x2){
	matrix<Type> ans(2,2);
	ans(0,0) = exp(logsigma_x1);
	ans(0,1) = 0;
	ans(1,0) = 0;
	ans(1,1) = exp(logsigma_x2);
	return ans;
}
template<class Type>
vector<Type> h(Type x1, Type x2){
	vector<Type> ans(2);
	ans(0) = x1;
	ans(1) = x2;
	return ans;
}
template<class Type>
matrix<Type> construct_H(matrix<Type> Xsp, Type a){
  int nn = 5;
  int m = 2;
  matrix<Type> H(m,nn);
  vector<Type> x0;
  for(int i=0;i<nn;i++){
    x0 = Xsp.col(i);
    H.col(i) = h(x0(0), x0(1));
  }
  return H;
}
template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(t);
	DATA_VECTOR(y1);
	DATA_VECTOR(y2);
	PARAMETER(loglambda);
	PARAMETER(logsigma_x1);
	PARAMETER(logsigma_x2);
	PARAMETER(logsigma_y1);
	PARAMETER(logsigma_y2);
	DATA_VECTOR(X0);
	DATA_MATRIX(P0);
	DATA_SCALAR(dt);
	/* DEFINE VARIABLES*/
	int __n = 2;
	int __nn = 5;
	int __m = 2;
	int __N,__s;
	Type __nll = 0;

	vector<Type> __y(__m);
	vector<Type> __NAs,__ynew,__e,__X1;
	matrix<Type> __Xsp0,__S0,__G, __F, __S0Inv, __M, __S0PhiM, __Frhs1(__n,__nn), __Frhs0, __Frhs;
	matrix<Type> __E, __H, __Syy, __SyyInv, __Sxy, __K, __P1;
	vector<vector<Type>> xPost(t.size());
	vector<matrix<Type>> pPost(t.size());
	vector<vector<Type>> OneStepErrors(t.size());
	vector<matrix<Type>> OneStepErrorsCovariance(t.size());
	xPost(0) = X0;
	pPost(0) = P0;

	matrix<Type> __V(__m,__m);
	__V.setZero();
	__V.diagonal() << exp(logsigma_y2), exp(logsigma_y1);

	/*CONSTRUCT WEIGHTS*/
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
  __W = (__I-__Wmm) * __Wm_diag * (__I-__Wmm).transpose();


	/* MAIN FOR-LOOP - FILTERING/UPDATE*/
	 __S0 = P0.llt().matrixL();
	 __Xsp0 = construct_Xsp(X0,__S0);
	 for(int i=0;i<t.size()-1;i++){
		 __N = CppAD::Integer((t(i+1)-t(i))/dt);
		 for(int j=0;j<__N;j++){
			 __G  = g(logsigma_x1, logsigma_x2);
			 __F  = construct_F(__Xsp0,loglambda);
			 __S0Inv = __S0.inverse();
			 __M = __S0Inv * (__Xsp0 * __W * __F.transpose() + __F * __W * __Xsp0.transpose() + __G*__G.transpose()) * __S0Inv.transpose();
			 __S0PhiM = __S0 * Phi(__M);
			 __Frhs1.block(0,1,__n,__n) = __S0PhiM;
			 __Frhs1.block(0,__n+1,__n,__n) = -__S0PhiM;
			 __Frhs0 = (__F*__Wm).replicate(1,__nn);
			 __Frhs = __Frhs0 + sqrt(3.0) * __Frhs1;
			 __Xsp0 += __Frhs * dt;
			 __S0 = ((__Xsp0 - __Xsp0.col(0).replicate(1,__nn))/sqrt(Type(3.0))).block(0,1,__n,__n);
			};
		__P1 = __S0 * __S0.transpose();
		__y << y1(i+1), y2(i+1);
		__NAs = isNA(__y);
		__s = CppAD::Integer(sum(__NAs));
		if( __s > 0 ){
			 __ynew = removeNAs(__s,__y,__NAs);
			 matrix<Type> __E0(__s,__m);
			 __E 	= constructE(__E0,__NAs);
			 __H = construct_H(__Xsp0,Type(0.0));
			 __e  = __ynew - __E * (__H * __Wm);
			 __Syy  = __E * (__H * __W * __H.transpose() + __V) * __E.transpose();
			 __Sxy  = __Xsp0 * __W * __H.transpose() * __E.transpose();
			 __SyyInv  = __Syy.inverse();
			 __K = __Sxy * __SyyInv;
			 __X1 = __Xsp0 * __Wm + __K * __e;
			 __P1 = __S0 * __S0.transpose() - __K * __Syy * __K.transpose();
			 __S0 = __P1.llt().matrixL();
			 __Xsp0 = construct_Xsp(__X1,__S0);
			 __nll += 0.5*atomic::logdet(__Syy) + 0.5*(__e*(__SyyInv*__e)).sum() + 0.5*log(2*M_PI)*asDouble(__s);
			 OneStepErrors(i) = __e;
			 OneStepErrorsCovariance(i) = __Syy;
		};
		 xPost(i+1) = __Xsp0.col(0);
		 pPost(i+1) = __P1;
	};
	REPORT(xPost);
	REPORT(pPost);
	REPORT(OneStepErrors);
	REPORT(OneStepErrorsCovariance);
	return __nll;
};
