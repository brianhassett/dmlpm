#' Compute the full log-posterior
#'
#' @param Y
#' @param Z
#' @param alpha
#' @param sigma2_0
#' @param sigma2_a
#' @param sigma2_z
#'
#' @return
#' @export
#'
#' @examples
log_posterior <- function(Y, Z, alpha, sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z) {

  # Data params:
  n <- nrow(Y)
  TT <- dim(Y)[3]
  K <- dim(Y)[4]

  # Compute pairwise distance matrix for each time point (3D array, dim = (n, n, T)):
  D <- array(NA, dim = c(n, n, TT))
  for(t in 1:TT) {
    D[ , , t] <- Dist(Z[ , , t])
  }

  # Log-likelihood:
  l <- 0
  for(k in 1:K) {
    for(t in 1:TT) {
      Y_tk <- Y[ , , t, k] # extract t_th slice of k_th layer
      # Dec 7 2024 - just add an index for t here for the dynamic alpha:
      l <- l + sum(upper_tri(Y_tk)*(alpha[k, t] - upper_tri(D[ , , t])) - exp(alpha[k, t] - upper_tri(D[ , , t])) - lfactorial(upper_tri(Y_tk)))
    }
  }

  # Log posterior:
  log_post <-
    l +                                                                # log-likelihood
    #K*log(1/sqrt(2*pi*sigma2_a)) - sum((alpha^2)/(2*sigma2_a)) +       # log-prior for alpha (Dec 7 2024: old for static alpha)

    # New prior for dynmaic alpha (this is based on the one for z_t, which is also a RW and very similar:
    K*log(1/sqrt(2*pi*sigma2_a)) - sum((alpha[ , 1]^2)/(2*sigma2_a)) + # log-prior for alpha @ t = 1 - The "initial mean" is zero for now
    K*(TT - 1)*log(1/sqrt(2*pi*sigma2_a)) - (1/(2*sigma2_a))*sum((alpha[ , 2:TT] - alpha[ , 1:(TT-1)])^2) + # log-prior for alpha @ t > 1

    n*log(1/(2*pi*sigma2_0)) - (1/(2*sigma2_0))*sum(Z[ , , 1]^2) +     # log-prior for Z_1 (t = 1)
    n*(TT - 1)*log(1/(2*pi*sigma2_z)) - (1/(2*sigma2_z))*sum((Z[ , , 2:TT] - Z[ , , 1:(TT-1)])^2) # log-prior for Z_t (t = 2 to t = T)

  return(log_post)
}


#' Compute the contribution to the log-posterior after updating at a given time t and layer k
#'
#' @param Y_tk
#' @param D_t
#' @param alpha
#' @param sigma2_0
#' @param sigma2_a
#' @param sigma2_z
#' @param k
#' @param t
#'
#' @return
#' @export
#'
#' @examples
log_posterior_alpha_kt <- function(Y_tk, D_t, alpha, sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z, k = 1, t = 1) {

  # Data params:
  n <- nrow(Y)
  TT <- dim(Y)[3]
  K <- dim(Y)[4]

  D_t_upper <- upper_tri(D_t)

  # Log-likelihood:
  l <- sum(Y_tk*(alpha[k, t] - D_t_upper) - exp(alpha[k, t] - D_t_upper))

  # Log posterior:
  log_post <-
    l -                                                      # log-likelihood
    (alpha[k, 1]^2)/(2*sigma2_a) - # log-prior for alpha @ t = 1 - The "initial mean" is zero for now
    (1/(2*sigma2_a))*sum((alpha[k, 2:TT] - alpha[k, 1:(TT-1)])^2)  # log-prior for alpha @ t > 1

  return(log_post)
}
