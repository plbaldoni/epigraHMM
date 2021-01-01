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

// E-step of HMM (forward-backward probability + posterior probability calculation)
//[[Rcpp::export]]
void saveMixtureProb(arma::mat eta,
                     Rcpp::StringVector hdf5){
    
    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    
    // Normalizing
    eta = normalise(eta,1,1);
    
    // Saving mixture probabilities to file
    eta.save(hdf5_name(vstrings[0], "mixtureProb",hdf5_opts::replace));
}