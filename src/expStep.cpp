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
void expStep(arma::mat counts,
             arma::vec pi,
             arma::mat gamma,
             arma::mat logf,
             Rcpp::StringVector nameForwardProb,
             Rcpp::StringVector nameBackwardProb,
             Rcpp::StringVector nameMarginalProb,
             Rcpp::StringVector nameJointProb){
    
    int M = counts.n_rows; //Number of observations
    int K = pi.n_elem; //Number of states
    
    std::vector<std::string> vstringsFP(1);
    std::vector<std::string> vstringsBP(1);
    std::vector<std::string> vstringsP1(1);
    std::vector<std::string> vstringsP2(1);
    
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
        // Auxiliary quantities for FB algorithm
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
    
    // Saving forward-backward probabilities
    vstringsFP[0] = nameForwardProb(0);
    vstringsBP[0] = nameBackwardProb(0);
    logFP.save(vstringsFP[0]);
    logBP.save(vstringsBP[0]);
    
    // Calculating posterior probabilities now
    arma::mat logProb1(M,K); //Window-based posterior probabilities
    arma::mat logProb2(M,K*K); //Window-based joint posterior probabilities

    // Avoiding underflow
    arma::mat d(1,K);
    d.row(0) = logFP.row(M-1);
    double maxd = max(d.row(0));
    
    // Auxiliary matrices for joint probabilities
    IntegerVector KK = seq_len(K*K)-1;
    int KKsize = KK.size();
    
    IntegerVector K_each = rep_each(seq_len(K)-1,K);
    IntegerVector K_rep = rep(seq_len(K)-1,K); 
    
    // Looping through windows
    for(int j = 0; j < M;j++){
        logProb1.row(j) = logFP.row(j) + logBP.row(j) - (maxd + log(sum(exp(d.row(0)-maxd))));
        for(int k = 0; k < KKsize; k++){
            if(j >=1){
                logProb2(j,k) = logFP(j-1,K_each(k)) + log(gamma(K_each(k),K_rep(k))) + logf(j,K_rep(k)) + logBP(j,K_rep(k)) - (maxd + log(sum(exp(d.row(0)-maxd))));
            }
        }
    }
    
    // Saving the posterior probabilities
    vstringsP1[0] = nameMarginalProb(0);
    logProb1.save(vstringsP1[0]);
    vstringsP2[0] = nameJointProb(0);
    logProb2.save(vstringsP2[0]);
}