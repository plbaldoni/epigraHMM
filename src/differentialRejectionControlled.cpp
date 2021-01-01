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
arma::field<arma::mat> differentialRejectionControlled(Rcpp::StringVector hdf5,
                                arma::vec f,
                                double p,
                                int N) { 
    
    arma::mat logProb1;
    arma::mat reWeights;
    arma::mat eta;

    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    logProb1.load(hdf5_name(vstrings[0], "logProb1"));
    eta.load(hdf5_name(vstrings[0], "mixtureProb"));

    int K = logProb1.n_cols; //Number of states
    int B = eta.n_cols;
    int KB = K + B - 1;
    arma::field<mat> aggReWeights(KB,1);

    // Consensus background
    reWeights = repmat(reweight(exp(logProb1.col(0)),p),N,1);
    aggReWeights(0,0) = aggregate(reWeights.col(0),f);

    // Differential
    for(int b = 0; b < B; b++){
        reWeights = repmat(reweight(exp(logProb1.col(1))%eta.col(b),p),N,1);
        aggReWeights(b+1,0) = aggregate(reWeights.col(0),f);
    }

    // Consensus enrichment
    reWeights = repmat(reweight(exp(logProb1.col(2)),p),N,1);
    aggReWeights(KB-1,0) = aggregate(reWeights.col(0),f);
    
    return aggReWeights;
}