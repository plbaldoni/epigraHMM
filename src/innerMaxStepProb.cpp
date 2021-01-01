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

// inner M-step (maximization w.r.t. mixture probabilities)
//[[Rcpp::export]]
arma::vec innerMaxStepProb(Rcpp::StringVector hdf5){
    
    arma::mat eta; // mixture probabilities
    arma::mat logProb1; // marginal probabilities

    // Loading probabilities
    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    eta.load(hdf5_name(vstrings[0], "mixtureProb"));
    logProb1.load(hdf5_name(vstrings[0], "logProb1"));

    // Number of mixture components
    int B = eta.n_cols;
    
    // MLE
    arma::vec delta(B);
    
    // Looping
    for(int b = 0; b < B; b++){
        delta(b) = accu(exp(logProb1.col(1))%eta.col(b))/accu(exp(logProb1.col(1)));
    }
    
    return delta;
}