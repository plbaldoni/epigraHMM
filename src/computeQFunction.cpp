//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Compute the Q-function of the EM-algorithm
//[[Rcpp::export]]
double computeQFunction(Rcpp::StringVector hdf5,
                           arma::vec pi,
                           arma::mat gamma) { 
    
    // Loading log-likelihood function, and posterior probabilities 
    arma::mat logf;
    arma::mat Prob1;
    arma::mat Prob2;
    
    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    logf.load(hdf5_name(vstrings[0], "logLikelihood"));
    Prob1.load(hdf5_name(vstrings[0], "logProb1"));
    Prob2.load(hdf5_name(vstrings[0], "logProb2"));
    
    Prob1 = exp(Prob1);
    Prob2 = exp(Prob2);
    
    double Q = sum(Prob1.row(0)*log(pi)) + accu(Prob1%logf) + accu(Prob2.rows(1,Prob2.n_rows - 1)*diagmat(log(vectorise(gamma,1))));
    
    return Q;
}