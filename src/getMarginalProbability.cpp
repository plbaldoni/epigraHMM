//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Get marginal probabilities saved in disk
//[[Rcpp::export]]
arma::mat getMarginalProbability(Rcpp::StringVector hdf5) { 
    
    arma::mat logProb1;
    
    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    logProb1.load(hdf5_name(vstrings[0], "logProb1"));
    
    return exp(logProb1);
}