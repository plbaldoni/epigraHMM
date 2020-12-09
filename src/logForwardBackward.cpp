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

//Forward-Backward Probabilities for MHMM
//[[Rcpp::export]]
void logForwardBackward(arma::mat counts, arma::vec pi, arma::mat gamma, arma::mat logf,
                        Rcpp::StringVector nameF, Rcpp::StringVector nameB){
    //logf = log-likelihood matrix
    //pi = probability vector
    //gamma = transition probability
    
    int M = counts.n_rows; //Number of observations
    int K = pi.n_elem; //Number of states
    std::vector<std::string> vstringsF(1);
    std::vector<std::string> vstringsB(1);
    
    // Creating likelihood
    // For a vectorized version, see https://stackoverflow.com/questions/29430726/vectorised-rcpp-random-binomial-draws
    
    arma::mat logFP(M,K); // log-forward probabilities
    arma::mat logBP(M,K); // log-backward probabilities
    arma::rowvec headRow(K);
    arma::rowvec tailRow(K);
    
    arma::mat a(1,K);
    double maxa = 0;
    arma::mat b(1,K);
    double maxb = 0;
    
    // Forward Backward algorithm
    
    logFP.row(0) = log(pi.t()) + logf.row(0);
    logBP.row(M-1).zeros();
    for(int j = 1; j < M; j++){
        headRow = logf.row(j);
        tailRow = logf.row(M-j);
        for(int k = 0; k < K; k++){
            // Forward algorithm
            a.row(0) = logFP.row(j-1) + log(trans(gamma.col(k))) + headRow(k);
            maxa = max(a.row(0));
            logFP(j,k) = maxa + log(sum(exp(a.row(0)-maxa)));
            
            // Backward algorithm
            b.row(0) = logBP.row(M-j) + log(gamma.row(k)) + tailRow;
            maxb = max(b.row(0));
            logBP(M-j-1,k) = maxb + log(sum(exp(b.row(0)-maxb)));
        }
    }
    
    //Saving logFB probabilities for i-th cell
    vstringsF[0] = nameF(0);
    vstringsB[0] = nameB(0);
    logFP.save(vstringsF[0]);
    logBP.save(vstringsB[0]);
}