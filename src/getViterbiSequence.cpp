//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
using namespace Rcpp;
using namespace arma;
using namespace H5;

// Compute the Viterbi sequence
//[[Rcpp::export]]
arma::vec getViterbiSequence(Rcpp::StringVector hdf5, arma::vec pi, arma::mat gamma) { 
    
    arma::mat logFP;
    
    std::vector<std::string> vstrings(1);
    vstrings[0] = hdf5(0);
    logFP.load(hdf5_name(vstrings[0], "logFP"));
    
    int M = logFP.n_rows;
    int K = logFP.n_cols;
    
    arma::mat LOGV(M,K);
    arma::mat V(M,K);
    
    LOGV.row(0) = log(pi.t()) + logFP.row(0);
    
    arma::mat COND(K,K);
    arma::vec MAXCOND(K);
    arma::vec S(M);
    arma::vec AUXS(K);
    
    for(int j = 0; j < (M-1); j++){
        for(int k = 0; k < K; k++){
            COND.row(k) = LOGV.row(j) + log(gamma.col(k).t());
            MAXCOND[k] = max(COND.row(k));
        }
        LOGV.row(j+1) = MAXCOND.t() + logFP.row(j+1);
    }

    S[M-1] = LOGV.row(M-1).index_max();

    for(int j = (M-2); j >= 0; j--){
        AUXS = LOGV.row(j).t() + log(gamma.col(S[j+1]));
        S[j] = AUXS.index_max();
    }

    return S;
}