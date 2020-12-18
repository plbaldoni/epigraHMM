//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
#include "reweight.h"
#include "aggregate.h"
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Applies reweight.cpp and aggregate.cpp for rejection-controlled posterior probabilities
//[[Rcpp::export]]
arma::mat rejectionControlled(NumericVector x,arma::vec f,double p) { 
    
    arma::vec reWeights = reweight(x,p); // Reweights posterior probabilities
    return aggregate(reWeights,f); // Aggregating
}