#include <RcppArmadillo.h>
using namespace Rcpp;

// For the Dynamic Multiview LPM Mrch 11th, see lpm_dynamic_multiview_3.Rmd for first use

// Dec 7 2024: added dynamic alpha stuff (fixed variance)


// Function to calculate pairwise distance matrix using RcppArmadillo
// [[Rcpp::export]]
arma::rowvec pairwise_distance_vector_cpp(arma::mat Z, int i) {
  int n = Z.n_rows;
  arma::rowvec d(n);

  for (int j = 0; j < n; j++) {
    if (j != i) {
      arma::rowvec diff = (Z.row(i) - Z.row(j));
      d(j) = sqrt(arma::dot(diff, diff));
    }
  }

  return d;
}

// Rcpp function to compute log posterior
// [[Rcpp::export]]
double log_posterior_it_cpp(const arma::field<arma::cube>& Y, arma::cube Z, const NumericMatrix& alpha, int i, int t, double sigma2_0 = 1, double sigma2_a = 1, double sigma2_z = 1) {

  // Log-likelihood for the i-th node:
  int T = Z.n_slices;
  int K = Y.n_elem;

  // Log prior for Z_i:
  double log_prior_z_i;
  if(t == 0) { // t = 1
    log_prior_z_i = (1 / (2 * sigma2_0)) * (pow(Z(i, 0, 0), 2) + pow(Z(i, 1, 0), 2)) +
      (1 / (2 * sigma2_z)) * (pow(Z(i, 0, 1) - Z(i, 0, 0), 2) + pow(Z(i, 1, 1) - Z(i, 1, 0), 2));
  } else if(t == (T - 1)) { // t = T
    log_prior_z_i = (1 / (2 * sigma2_z)) * (pow(Z(i, 0, t) - Z(i, 0, t - 1), 2) + pow(Z(i, 1, t) - Z(i, 1, t - 1), 2));
  } else { // 1 < t < T
    log_prior_z_i = (1 / (2 * sigma2_z)) * (pow(Z(i, 0, t) - Z(i, 0, t - 1), 2) + pow(Z(i, 1, t) - Z(i, 1, t - 1), 2)) +
      (1 / (2 * sigma2_z)) * (pow(Z(i, 0, t + 1) - Z(i, 0, t), 2) + pow(Z(i, 1, t + 1) - Z(i, 1, t), 2));
  }

  // Likelihood:
  arma::rowvec d = pairwise_distance_vector_cpp(Z.slice(t), i);
  d.shed_col(i); // drop the i-th entry, equivalent to d[-i]

  double l = 0;
  for (int k = 0; k < K; k++) {

    arma::rowvec Y_row_i = Y(k).slice(t).row(i);
    Y_row_i.shed_col(i);

    l += sum(Y_row_i % (alpha(k, t) - d) - exp(alpha(k, t) - d));
  }

  // Log posterior:
  double log_post_i = l - log_prior_z_i;

  return log_post_i;
}
