
#' Print Procrustes correlation of two sets of latent positions
#'
#' @param Z1 # n x p x T array
#' @param Z2 # n x p x T array
#'
#' @return
#' @export
#'
#' @examples
procrustes_cor <- function(Z1, Z2) {
  TT <- dim(Z1)[3]
  walk(1:TT, function(t) print(paste0("t = ", t, ": ", round(vegan::protest(Z1[ , , t], Z2[ , , t])$t0, 3))))
}

## TODO
## Fix this so that fit$Z_trans is passed, rather than using it in the function.

procrustes_cor_stacked <- function() {
  # Procrustes correlation of the stacked matrices:
  # (first I need to take the posterior mean of the stacked and transformed samples, the Z_trans matrix was already made above)
  Z_trans_mean <- matrix(NA, nrow = dim(fit$Z_trans)[1], ncol = dim(fit$Z_trans)[2])
  for(i in 1:(n*TT)) {
    for(j in 1:2) {
      Z_trans_mean[i, j] <- mean(fit$Z_trans[i, j, ])
    }
  }

  # Construct stacked matrix of the true Z matrices:
  Z_stack_true <- matrix(NA, nrow = n*TT, ncol = d)
  for(t in 1:TT) {
    Z_stack_true[((t - 1)*n + 1):(t*n), ] <- Z[ , , t] # These indexes look weird but they are correct:
  }

  # Compare to the stacked true matrices:
  paste0("Correlation of stacked matrices (nT x 2): ", round(vegan::protest(Z_stack_true, Z_trans_mean)$t0, 3))
}
