# Function adapted from latentnet to compute geodesic distance matric:
#' Title
#'
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
ergmm.geodesicmatrix <- function(Y) {
  # For the purpose of geodesic distance, dichotomize the network about its mean (Brian: weighted networks only)
  #Y <- Y > mean(Y, na.rm=TRUE)
  Y[is.na(Y)] <- 0
  Y <- Y | t(Y) # Brian: double check what this does
  sna::geodist(Y, count.paths = FALSE, inf.replace = nrow(Y))$gdist
}

# Function to perform MDS adapted from latentnet:
#' Title
#'
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
mds_initialise <- function(Y) {

  # Compute distance natrix:
  D <- ergmm.geodesicmatrix(Y)
  D[is.infinite(D)] <- 2*n

  # Apply MDS to obtain initial postions:
  Z_0 <- cmdscale(D, 2)

  # Checks:
  i.keep <- mahalanobis(Z_0, 0, cov(Z_0)) < 20
  Z_0[!i.keep,] <- 0

  # Return initial positions:
  return(Z_0)
}
