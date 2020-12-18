#ifndef REWEIGHT_H
#define REWEIGHT_H

#include <RcppArmadillo.h>
#include <numeric>
#include <Rcpp.h>
#include <H5Cpp.h>
#include "rbinomVectorized.h"
arma::vec reweight(arma::vec x,double p);

#endif