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
arma::field<arma::mat> rejectionControlled(Rcpp::StringVector nameMarginalProb,arma::vec f,double p) { 
    arma::mat logProb1;
    arma::vec reWeights;

    std::vector<std::string> vstrings(1);
    vstrings[0] = nameMarginalProb(0);
    logProb1.load(vstrings[0]);
    
    int K = logProb1.n_cols; //Number of states
    arma::field<mat> aggReWeights(K,1);

    for(int k = 0; k < K; k++){
        reWeights = reweight(exp(logProb1.col(k)),p);
        aggReWeights(k,0) = aggregate(reWeights,f);
    }

    return aggReWeights;
}