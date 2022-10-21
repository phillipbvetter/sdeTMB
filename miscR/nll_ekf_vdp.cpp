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
vector<Type> f(Type loglambda, Type x1, Type x2){
	vector<Type> ans(2);
	ans(0) = x2;
	ans(1) = x2 * (1 - pow(x1, 2)) * exp(loglambda) - x1;
	return ans;
}
template<class Type>
matrix<Type> dfdx(Type loglambda, Type x1, Type x2){
	matrix<Type> ans(2,2);
	ans(0,0) = 0;
	ans(0,1) = 1;
	ans(1,0) = -(1 + 2 * (x1 * x2 * exp(loglambda)));
	ans(1,1) = (1 - pow(x1, 2)) * exp(loglambda);
	return ans;
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
matrix<Type> dhdx(Type a){
	matrix<Type> ans(2,2);
	ans(0,0) = 1;
	ans(0,1) = 0;
	ans(1,0) = 0;
	ans(1,1) = 1;
	return ans;
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
	int __n = 2;
	int __m = 2;
	int __k = CppAD::Integer((diff(t)/dt).sum() + 1);
	int __N;
	int __s;
	Type __nll = 0;

	vector<Type> __x0 = X0;
	matrix<Type> __p0 = P0;
	vector<Type> __x1 = __x0;
	matrix<Type> __p1 = __p0;
	vector<vector<Type>> xPrior(t.size());
	vector<vector<Type>> xPost(t.size());
	vector<vector<Type>> xPriorPost(2*t.size());
	vector<matrix<Type>> pPrior(t.size());
	vector<matrix<Type>> pPost(t.size());
	vector<matrix<Type>> pPriorPost(2*t.size());
	vector<vector<Type>> xPrior_all(__k);
	vector<matrix<Type>> pPrior_all(__k);
	vector<vector<Type>> F_all(__k);
	vector<matrix<Type>> A_all(__k);
	vector<matrix<Type>> G_all(__k);
	vector<vector<Type>> OneStepErrors(t.size());
	vector<matrix<Type>> OneStepErrorsCovariance(t.size());
	xPrior(0) = __x0;
	pPrior(0) = __p0;
	xPost(0) = __x0;
	pPost(0) = __p0;
	xPriorPost(0) = __x0;
	xPriorPost(1) = __x0;
	pPriorPost(0) = __p0;
	pPriorPost(1) = __p0;
	matrix<Type> __C,__R,__K,__E,__C2,__V2,__Ri,__A,__G;
	vector<Type> __e0,__e,__y(__m),__NAs,__F,__ynew,__h;
	matrix<Type> __I(__n,__n);
	__I.setIdentity();
	matrix<Type> __V(__m,__m);
	__V.setZero();
	__V.diagonal() << exp(logsigma_y2), exp(logsigma_y1);

	for(int i=0;i<t.size()-1;i++){
		__N = CppAD::Integer((t(i+1)-t(i))/dt);
		for(int j=0;j<__N;j++){
			__A  = dfdx(loglambda, __x0(0), __x0(1));
			__F  = f(loglambda, __x0(0), __x0(1));
			__G  = g(logsigma_x1, logsigma_x2);
			__x1 = __x0 + __F * dt;
			__p1 = __p0 + ( __A*__p0 + __p0*__A.transpose() + __G*__G.transpose() ) * dt;
			__x0 = __x1;
			__p0 = __p1;
			xPrior_all(i*__N+j) = __x0;
			pPrior_all(i*__N+j) = __p0;
			F_all(i*__N+j) = __F;
			A_all(i*__N+j) = __A;
			G_all(i*__N+j) = __G;
		}
		xPrior(i+1) = __x0;
		pPrior(i+1) = __p0;
		xPriorPost(2*i+2) = __x0;
		pPriorPost(2*i+2) = __p0;
		__y << y1(i+1), y2(i+1);
		__NAs = isNA(__y);
		__s = CppAD::Integer(sum(__NAs));
		if( __s > 0 ){
			__ynew = removeNAs(__s,__y,__NAs);
			matrix<Type> __E0(__s,__m);
			__E 	= constructE(__E0,__NAs);
			__h  = h(__x0(0), __x0(1));
			__e  = __ynew - __E*__h;
			__C  	= dhdx(Type(0.0));
			__C2 	= __E*__C;
			__V2  = __E*__V*__E.transpose();
			__R 	= __C2*__p0*__C2.transpose() + __V2;
			__Ri  = __R.inverse();
			__K 	= __p0 * __C2.transpose() * __Ri;
			__x0  = __x0 + __K*__e;
			__p0  = (__I-__K*__C2)*__p0*(__I-__K*__C2).transpose() + __K*__V2*__K.transpose();
			__nll += 0.5*atomic::logdet(__R) + 0.5*(__e*(__Ri*__e)).sum() + 0.5*log(2*M_PI)*asDouble(__s);
			OneStepErrors(i+1) = __e;
			OneStepErrorsCovariance(i+1) = __R;
		}
		xPost(i+1) = __x0;
		pPost(i+1) = __p0;
		xPriorPost(2*i+3) = __x0;
		pPriorPost(2*i+3) = __p0;
	}

REPORT(OneStepErrors);
REPORT(OneStepErrorsCovariance);
REPORT(xPrior);
REPORT(xPost);
REPORT(xPriorPost);
REPORT(xPrior_all);
REPORT(pPrior);
REPORT(pPost);
REPORT(pPriorPost);
REPORT(pPrior_all);
REPORT(F_all);
REPORT(A_all);
REPORT(G_all);
return __nll;
}
