//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
#include "rbinomVectorized.h"
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Enable lambda expressions.... 
// [[Rcpp::plugins(cpp11)]]

// Re-weights using a rejection-controlled approach
//[[Rcpp::export]]
arma::vec reweight(arma::vec x,double p) { 
    
    arma::uvec ids = find(x < p); // Find indices
    x.elem(ids) = rbinomVectorized(x.elem(ids)/p)*p; // Assign value to condition 
    
    return(x);
}