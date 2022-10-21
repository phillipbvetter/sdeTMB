#include <TMB.hpp>
using namespace density;
template<class Type>
bool isNA(Type x){
	return R_IsNA(asDouble(x));
}
template<class Type>
matrix<Type> VarMat(Type logsigma_x, Type logtheta){
	matrix<Type> G(2,2);
	G(0,0) = exp(logtheta);
	G(0,1) = 0 + exp(logsigma_x) * exp(logsigma_x);
	G(1,0) = 0;
	G(1,1) = -exp(logtheta);
	return G;
}
template<class Type>
matrix<Type> MeanMat(Type logtheta, Type mu){
	matrix<Type> A(2,2);
	A(0,0) = -exp(logtheta);
	A(0,1) = mu * exp(logtheta);
	A(1,0) = 0;
	A(1,1) = 0;
	return A;
}
template<class Type>
Type objective_function<Type>::operator() ()
{
	DATA_VECTOR(t);
	DATA_VECTOR(y);
	PARAMETER_VECTOR(x);
	PARAMETER(logtheta);
	PARAMETER(mu);
	PARAMETER(logsigma_x);
	PARAMETER(logsigma_y);
	Type __nll = 0;
	vector<Type> __Xi(1);
	vector<Type> __Xip1(1);
	vector<Type> __Z(1);
	matrix<Type> __MeanAug0(2,2);
	matrix<Type> __VarAug0(2,2);
	matrix<Type> __MeanAug(2,2);
	matrix<Type> __VarAug(2,2);
	matrix<Type> __Ahat(1,1);
	vector<Type> __Bhat(1);
	matrix<Type> __Qhat(1,1);
	Type __dt;
	__dt = t(1)-t(0);
	__MeanAug0 = MeanMat(logtheta, mu)*__dt;
	__VarAug0 = VarMat(logsigma_x, logtheta)*__dt;
	__MeanAug = expm(__MeanAug0);
	__VarAug = expm(__VarAug0);
	__Ahat = __MeanAug.block(0,0,1,1);
	__Bhat = __MeanAug.col(1).head(1);
	__Qhat = __VarAug.block(1,1,1,1).transpose() * __VarAug.block(0,1,1,1);

	for(int i=0;i<t.size()-1;i++){
		__Xi << x(i);
		__Xip1 << x(i+1);
		__Z = __Xip1 - (__Ahat * __Xi + __Bhat);
		__nll += MVNORM(__Qhat)(__Z);
	}
	for(int i=0;i<y.size(); ++i){
		if(!isNA(y(i))){
			__nll -= dnorm(y(i),x(i),sqrt(exp(logsigma_y)),true);
		}
	}
return(__nll);
}
