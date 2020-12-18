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

int increment_maybe(int value, double vec_i){
    return vec_i == 0 ? value : ( value +1 ) ;  
}

// Aggregate a numeric vector by group
// From https://stackoverflow.com/questions/16975034/rcpp-equivalent-for-rowsum
//[[Rcpp::export]]
arma::mat aggregate(arma::vec x,arma::vec f) { 
    std::vector<double> vec(10) ;
    vec.reserve(1000); 
    int n=x.size(); 
    for( int i=0; i<n; i++){
        int index=f[i]; 
        while( index >= vec.size() ){
            vec.resize( vec.size() * 2 ) ;    
        }
        vec[ index ] += x[i] ;
    }
    // count the number of non zeros
    int s = std::accumulate( vec.begin(), vec.end(), 0, increment_maybe) ; 
    arma::vec result(s) ;
    arma::vec names(s) ;
    
    std::vector<double>::iterator it = vec.begin() ;
    for( int i=0, j=0 ; j<s; j++ ,++it, ++i ){
        // move until the next non zero value
        while( ! *it ){ i++ ; ++it ;}
        result[j] = *it ;
        names[j]  = i ;
    }

    return arma::join_horiz(result,names);
}