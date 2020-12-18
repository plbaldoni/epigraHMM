#ifndef AGGREGATE_H
#define AGGREGATE_H

#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
using namespace Rcpp;
using namespace arma;
using namespace H5;
arma::mat aggregate(arma::vec x,arma::vec f);

#endif