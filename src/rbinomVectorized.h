#ifndef RBINOMVECTORIZED_H
#define RBINOMVECTORIZED_H

#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
arma::vec rbinomVectorized(arma::vec prob);

#endif