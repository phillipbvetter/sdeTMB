#include <TMB.hpp>
using namespace density;
template<class Type>
matrix<Type> VarMat(Type logsigmaX, Type logtheta){
	matrix<Type> G(2,2);
	G(0,0) = exp(logtheta);
	G(0,1) = exp(logsigmaX)*exp(logsigmaX);
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
	DATA_VECTOR(Y);
	DATA_VECTOR(iobsY);
	PARAMETER_VECTOR(X);
	PARAMETER(logtheta);
	PARAMETER(mu);
	PARAMETER(logsigmaX);
	PARAMETER(logsigmaY);
	Type nll = 0;
	vector<Type> Xi(1);
	vector<Type> Xip1(1);
	vector<Type> Z(1);
	matrix<Type> MeanAug0(2,2);
	matrix<Type> VarAug0(2,2);
	matrix<Type> MeanAug(2,2);
	matrix<Type> VarAug(2,2);
	matrix<Type> Ahat(1,1);
	vector<Type> Bhat(1);
	matrix<Type> Qhat(1,1);
	Type dt;
