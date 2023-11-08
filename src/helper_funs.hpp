#include <Rcpp.h>
#include <RcppEigen.h>
#include <Ziggurat.h>
using namespace Rcpp;
using namespace Eigen;
static Ziggurat::Ziggurat::Ziggurat zigg; //zigg.norm() draws from a normal distribution?

//#ifndef _HELPERFUNS_ 
//#define _HELPERFUNS_

//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
double invlogit(double x);

//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
Eigen::MatrixXd construct_permutation_matrix(int number_of_available_obs, int number_of_obs_eqs, Eigen::VectorXi bool_is_not_na_obsVec);

//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
Eigen::VectorXd remove_NAs(Eigen::VectorXd obsVec, int number_of_available_obs, Eigen::VectorXi bool_is_not_na_obsVec);

//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
template <typename T1, typename T2, typename T3>
List ode_integrator(
  T1 f__, 
  T2 g__,
  T3 dfdx__, 
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec, 
  Eigen::VectorXd parVec, 
  Eigen::VectorXd inputVec, 
  double timestep, 
  int ode_solver)
{

  // Call drift and convert from Rcpp::NumericMatrix to Eigen::MatrixXd
  Rcpp::NumericVector F = f__(stateVec, parVec, inputVec);
  Eigen::Map<Eigen::VectorXd> F_eigen = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(F);
  Eigen::MatrixXd P_rhs = cov_ode_1step(g__, dfdx__, covMat, stateVec, parVec, inputVec);


  // Initialize return values
  Eigen::VectorXd X1;
  X1.setZero();
  Eigen::MatrixXd P1;
  P1.setZero();

  // Forward Euler
  if(ode_solver == 1){
    X1 = stateVec + F_eigen * timestep;
    P1 = covMat + P_rhs * timestep;
  }
  if(ode_solver==2){}

  // Return
  return Rcpp::List::create(Named("X1")=X1, Named("P1")=P1);
}




//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
template<typename T2, typename T3>
Eigen::MatrixXd cov_ode_1step(
  T2 g__,
  T3 dfdx__,
  Eigen::MatrixXd covMat, 
  Eigen::VectorXd stateVec,
  Eigen::VectorXd parVec,
  Eigen::VectorXd inputVec)
{

  // Compute drift jacobian and diffusion 
  Rcpp::NumericMatrix A = dfdx__(stateVec, parVec, inputVec);
  Rcpp::NumericMatrix G = g__(stateVec, parVec, inputVec);

  // Convert from NumericMatrix to Eigen matrices
  Eigen::Map<Eigen::MatrixXd> A_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(A);
  Eigen::Map<Eigen::MatrixXd> G_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(G);

  // Calculate the right-hand side of covariance moment differential equation
  Eigen::MatrixXd cov_ode_1step = A_eigen * covMat + covMat * A_eigen.transpose() + G_eigen * G_eigen.transpose();

  // Return
  return cov_ode_1step;
}


//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
//////////////////// New Fun ////////////////////
template<typename T1, typename T2>
MatrixXd euler_maruyama_simulation(
  T1 f__, 
  T2 g__, 
  Eigen::MatrixXd stateMat, 
  Eigen::VectorXd parVec, 
  Eigen::VectorXd inputVec, 
  double timestep, 
  int nsims,
  int n,
  int m)
{
  // Returns a matrix with nsims rows and columns equal to the number of system states - each row corresponds to taking a single
  // Euler-Maruyama step
  MatrixXd stateMat_next(nsims, n);
  VectorXd stateVec, dW(m);
  NumericVector F;
  NumericMatrix G;
  double sqrt_timestep = sqrt(timestep);

  // Perform one-step simulation for each row in stateMat
  for(int i=0; i < nsims; i++){

    // extract state values
    stateVec = stateMat.row(i);

    // Convert F and G from Rcpp::Numeric to Eigen
    F = f__(stateVec, parVec, inputVec);
    Eigen::Map<Eigen::VectorXd> F_eigen = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(F);
    G = g__(stateVec, parVec, inputVec);
    Eigen::Map<Eigen::MatrixXd> G_eigen = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(G);

    // Generate dW vector by sampling from standard normal
    for(int i=0; i < m; i++){
      dW(i) = zigg.norm();
    }

    // Perform euler-maruyama step
    stateMat_next.row(i) =  stateVec + F_eigen * timestep + G_eigen * dW * sqrt_timestep;
  }

  // Return
  return stateMat_next;
}



//#endif // _HELPERFUNS_