#include <Rcpp.h>
#include <Ziggurat.h>
static Ziggurat::Ziggurat::Ziggurat zigg;
// [[Rcpp::export]]
Rcpp::NumericVector zrnorm(int n) {
Rcpp::NumericVector x(n); for (int i=0; i<n; i++) {
        x[i] = zigg.norm();
    }
return x; }
// [[Rcpp::export]]
void zsetseed(unsigned long int s) {
    zigg.setSeed(s);
return; }