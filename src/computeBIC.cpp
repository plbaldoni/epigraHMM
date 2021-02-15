//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Compute BIC for HMMs
//[[Rcpp::export]]
double computeBIC(Rcpp::StringVector hdf5,
                  int numPar,
                  int numSamples) { 
    
    // Loading log-likelihood function
    arma::mat logf;
    
    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    logf.load(hdf5_name(vstrings[0], "logFP"));
    
    int M = logf.n_rows;
    double C = max(logf.row(M-1));
    
    double BIC = -2*(C+log(sum(exp(logf.row(M-1)-C)))) + numPar*log(numSamples*M);
    
    return BIC;
}