#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double c_initialize1(NumericVector Data, IntegerVector DIMS, NumericVector Yy,
                     double XSCALE, double BETAIN, double BETAOUT, NumericVector WW) {

  arma::cube X(Data.begin(), DIMS[1], DIMS[2], DIMS[0]);
  arma::cube Y(Yy.begin(), DIMS[0], DIMS[0], DIMS[2]);

  double ret = 0, dx = 0, eta = 0;

  for (int tt = 0; tt < DIMS(2); tt++) {
    for (int i = 0; i < DIMS(0); i++) {
      for (int j = 0; j < DIMS(0); j++) {
        if (i != j) {
          dx = XSCALE * arma::norm(X.slice(i).col(tt) - X.slice(j).col(tt), 2);
          eta = (BETAIN * (1 - dx / WW[j]) + BETAOUT * (1 - dx / WW[i]));
          ret += Y.slice(tt)(i, j) * eta - log(1 + exp(eta));
        }
      }
    }
  }

  return ret;
}


// [[Rcpp::export]]
double c_initialize1_old(arma::cube X, arma::cube Y, double XSCALE, double BETAIN, double BETAOUT, arma::colvec WW) {
  int n = Y.n_rows;
  int T = Y.n_slices;
  double ret = 0, dx = 0, eta = 0;

  for(int tt = 0; tt < T; tt++) {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(i != j) {
          dx = XSCALE * arma::norm(X.slice(i).col(tt) - X.slice(j).col(tt), 2);
          eta = BETAIN * (1 - dx / WW(j)) + BETAOUT * (1 - dx / WW(i));
          ret += Y.slice(tt)(i, j) * eta - log(1 + exp(eta));
        }
      }
    }
  }
  return ret;
}


// [[Rcpp::export]]
arma::colvec c_initialize1_grad(NumericVector Data, IntegerVector DIMS, NumericVector Yy,
                                double XSCALE, double BETAIN, double BETAOUT, NumericVector WW) {

  arma::cube X(Data.begin(), DIMS[1], DIMS[2], DIMS[0]);
  arma::cube Y(Yy.begin(), DIMS[0], DIMS[0], DIMS[2]);

  int n = X.n_slices;  // Number of nodes
  int T = X.n_cols;    // Time steps

  arma::colvec ret(3, arma::fill::zeros); // Vector to store results

  for (int tt = 0; tt < T; tt++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j) {
          double dx = arma::norm(X.slice(i).col(tt) - X.slice(j).col(tt), 2);
          double eta = BETAIN * (1 - XSCALE * dx / WW(j)) + BETAOUT * (1 - XSCALE * dx / WW(i));
          double exp_eta = 1 / (1 + std::exp(-eta));

          ret(0) += dx * (BETAIN / WW(j) + BETAOUT / WW(i)) * (exp_eta - Y(i, j, tt));
          ret(1) += (1 - XSCALE * dx / WW(j)) * (Y(i, j, tt) - exp_eta);
          ret(2) += (1 - XSCALE * dx / WW(i)) * (Y(i, j, tt) - exp_eta);
        }
      }
    }
  }

  return ret;
}


// [[Rcpp::export]]
arma::vec c_initialize1_grad_old(arma::cube X, arma::cube Y, double XSCALE, double BETAIN, double BETAOUT, arma::colvec WW) {
  int n = Y.n_rows;
  int T = Y.n_slices;
  arma::vec ret(3, arma::fill::zeros);
  double dx = 0, eta = 0;

  for(int tt = 0; tt < T; tt++) {
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(i != j) {
          dx = arma::norm(X.slice(i).col(tt) - X.slice(j).col(tt), 2);
          eta = BETAIN * (1 - XSCALE * dx / WW(j)) + BETAOUT * (1 - XSCALE * dx / WW(i));

          double sigmoid = 1 / (1 + exp(-eta));
          ret(0) += dx * (BETAIN / WW(j) + BETAOUT / WW(i)) * (sigmoid - Y(i, j, tt));
          ret(1) += (1 - XSCALE * dx / WW(j)) * (Y(i, j, tt) - sigmoid);
          ret(2) += (1 - XSCALE * dx / WW(i)) * (Y(i, j, tt) - sigmoid);
        }
      }
    }
  }
  return ret;
}
