#include <TMB.hpp>
using namespace density;
template<class Type>
bool isNA(Type x){
	return R_IsNA(asDouble(x));
}
template<class Type>
matrix<Type> Gmat(Type logsigma_x1, Type logsigma_x2){
	matrix<Type> G(2,2);
	G(0,0) = exp(logsigma_x1);
	G(0,1) = 0;
	G(1,0) = 0;
	G(1,1) = exp(logsigma_x2);
	return(G);
}
template<class Type>
vector<Type> fvec(Type loglambda, Type x1, Type x2){
	vector<Type> f(2);
	f(0) = x2;
	f(1) = x2 * (1 - pow(x1, 2)) * exp(loglambda) - x1;
	return(f);
}
template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(t);
	DATA_VECTOR(y1);
	DATA_VECTOR(y2);
	PARAMETER_VECTOR(x1);
	PARAMETER_VECTOR(x2);
	PARAMETER(loglambda);
	PARAMETER(logsigma_x1);
	PARAMETER(logsigma_x2);
	PARAMETER(logsigma_y1);
	PARAMETER(logsigma_y2);
	int __n = 2;
	int __ng = 2;
	Type __nll = 0;
	vector<Type> __f(__n);
	vector<Type> __Xi(__n);
	vector<Type> __Xip1(__n);
	vector<Type> __Z(__n);
	matrix<Type> __G(__n,__ng);
	matrix<Type> __V(__n,__n);
	matrix<Type> __I(__n,__n);
	__I.setIdentity();
	__I *= 1e-8;
	Type __dt;
	for(int i=0;i<t.size()-1;i++){
		__dt = t(i+1) - t(i);
		__f = fvec(loglambda, x1(i), x2(i));
		__G = Gmat(logsigma_x1, logsigma_x2);
		__Xi << x1(i), x2(i);
		__Xip1 << x1(i+1), x2(i+1);
		__Z = __Xip1 - __Xi - __f * __dt;
		__V = (__G*__G.transpose() + __I) * __dt;
		__nll += MVNORM(__V)(__Z);
	}
	for(int i=0;i<y1.size(); ++i){
		if(!isNA(y1(i))){
			__nll -= dnorm(y1(i),x1(i),sqrt(exp(logsigma_y2)),true);
		}
	}
	for(int i=0;i<y2.size(); ++i){
		if(!isNA(y2(i))){
			__nll -= dnorm(y2(i),x2(i),sqrt(exp(logsigma_y1)),true);
		}
	}
return(__nll);
}
