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
//' Simulates a Markov Chain of length 'n' given a matrix of transition probabilities P
//' @param P a matrix of transition probabilities (row sums should be 1)
//' @param n an integer specifying thhe length of the simulated sequence
//'
//' @examples
//' #Example
//' simulateMarkovChain(matrix(c(0.99,0.01,0.01,0.99),2,2),100)
//' @export
//[[Rcpp::export]]
NumericVector simulateMarkovChain(NumericMatrix P,int n){
    NumericVector Z(n);
    Z[0] = 1;
    for(int i = 1; i < n; i ++)
    {
        double rand = R::runif(0,1);
        double sump = 0;
        int flag = 0;
        int k = 0;
        while(flag == 0){
            sump = sump + P((Z[i-1]-1),k);
            if(rand<=sump){
                Z[i] = (k+1);
                flag++;
            }
            k++;
        }
    }
    return Z;
}