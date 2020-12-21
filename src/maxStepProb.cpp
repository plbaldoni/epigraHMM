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

// M-step (maximization w.r.t. initial and transition probabilities)
//[[Rcpp::export]]
Rcpp::List maxStepProb(Rcpp::StringVector nameMarginalProb,
                       Rcpp::StringVector nameJointProb){
    
    arma::mat logProb1; // marginal probabilities
    arma::mat logProb2; // joint probabilities
    
    // Loading probabilities
    std::vector<std::string> vstringsProb1(1);
    vstringsProb1[0] = nameMarginalProb(0);
    logProb1.load(vstringsProb1[0]);
    
    std::vector<std::string> vstringsProb2(1);
    vstringsProb2[0] = nameJointProb(0);
    logProb2.load(vstringsProb2[0]);
    
    // Number of states & windows
    int M = logProb1.n_rows;
    int K = logProb1.n_cols;
    
    // Transition probability
    arma::mat gamma(K,K);
    double denom = 0;
    int idx1 = 0;
    int idx2 = 0;
    
    for(int k1 = 0; k1 < K; k1++){
        // Calculating the normalizing factor
        denom = accu(exp(logProb2.submat(1,idx1,M-1,idx1+K-1)));
        // Calculating MLE of transition probs.
        for(int k2 = 0; k2 < K; k2++){
            gamma(k1,k2) = accu(exp(logProb2.submat(1,idx2,M-1,idx2)));
            idx2 = idx2 + 1;
        }
        gamma.row(k1) = gamma.row(k1)/denom;
        idx1 = idx1 + K;
    }
    
    // Initial probability
    arma::vec pi = exp(logProb1.row(0).t());
    
    // Adding constraints to avoid underflow
    pi.elem(find(pi<0.000001)).fill(0.000001);
    pi = pi/sum(pi);
    gamma.elem(find(gamma<0.000001)).fill(0.000001);
    gamma = normalise(gamma,1,1);
    
    return Rcpp::List::create(Rcpp::Named("pi") = pi,
                              Rcpp::Named("gamma") = gamma);
}