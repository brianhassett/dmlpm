# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

c_initialize1 <- function(Data, DIMS, Yy, XSCALE, BETAIN, BETAOUT, WW) {
    .Call(`_dmlpm_c_initialize1`, Data, DIMS, Yy, XSCALE, BETAIN, BETAOUT, WW)
}

c_initialize1_old <- function(X, Y, XSCALE, BETAIN, BETAOUT, WW) {
    .Call(`_dmlpm_c_initialize1_old`, X, Y, XSCALE, BETAIN, BETAOUT, WW)
}

c_initialize1_grad <- function(Data, DIMS, Yy, XSCALE, BETAIN, BETAOUT, WW) {
    .Call(`_dmlpm_c_initialize1_grad`, Data, DIMS, Yy, XSCALE, BETAIN, BETAOUT, WW)
}

c_initialize1_grad_old <- function(X, Y, XSCALE, BETAIN, BETAOUT, WW) {
    .Call(`_dmlpm_c_initialize1_grad_old`, X, Y, XSCALE, BETAIN, BETAOUT, WW)
}

pairwise_distance_vector_cpp <- function(Z, i) {
    .Call(`_dmlpm_pairwise_distance_vector_cpp`, Z, i)
}

log_posterior_it_cpp <- function(Y, Z, alpha, i, t, sigma2_0 = 1, sigma2_a = 1, sigma2_z = 1) {
    .Call(`_dmlpm_log_posterior_it_cpp`, Y, Z, alpha, i, t, sigma2_0, sigma2_a, sigma2_z)
}

