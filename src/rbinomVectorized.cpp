//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Enable lambda expressions.... 
// [[Rcpp::plugins(cpp11)]]

// Vectorized version of rbinom
//[[Rcpp::export]]
arma::vec rbinomVectorized(arma::vec prob) { 
    arma::vec v(prob.size());
    std::transform(prob.begin(),prob.end(),v.begin(), [=](double p){ return R::rbinom(1,p);}); 
    return(v);
}