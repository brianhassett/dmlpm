#' Simulate networks from Poisson LPM given latent positions and intercept
#'
#' @param Z Matrix of latent positions
#' @param alpha Intercept parameter
#'
#' @return
#' @export
#'
#' @examples
LPM_gen_poisson <- function(Z, alpha) {
  N <- nrow(Z)
  adj <- matrix(0, N, N) # allocate an empty adjacency matrix
  for (i in 1:N) for (j in i:N) if (i != j)
  {
    d <- sqrt( sum( (Z[i , ] - Z[j, ])^2 ) ) # calculate the Euclidean distance between the two nodes
    lambda <- exp(alpha - d) # calculate rate
    adj[i,j] <- rpois(n = 1, lambda)
    adj[j,i] <- adj[i,j]
  }
  adj # return the adjacency matrix
}


