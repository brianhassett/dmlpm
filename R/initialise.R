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

#' Initialise Z based on Sewell and S&M (GMDS)
#'
#' @param Y
#' @param d
#'
#' @return
#' @export
#'
#' @examples
initialise_Z <- function(Y, d) { # d = dim of latent space, p in Sewell code

  # Brian
  n <- dim(Y)[1]
  TT <- dim(Y)[3]
  N <- 100 # Brian: N is iters, just set to 100 or anything I want since I only need the initialisations anyway

  w <- matrix(0,n,N)
  ###
  ###Weights (Brian: not sure what these are about, need to read on GMDS in S&M)
  ###
  for(tt in 1:TT){
    w[,1] <- w[,1] + apply(Y[,,tt,1],1,sum) +
      apply(Y[,,tt, 1],2,sum)
  }
  w[,1] <- w[,1]/sum(Y)/2
  if(sum(w==0)>0){
    w[,1] <- w[,1]+1e-5
    w[,1] <- w[,1]/sum(w[,1])
  }
  w[,1] <- w[,1]/sum(w[,1])

  ###
  ###Initial Latent Positions (GMDS, Sarkar and Moore, 2005)
  ###
  #p <- d # now an argument
  dissim <- array(0,dim=dim(Y[ , , , 1]))
  #for(tt in 1:TT) dissim[,,tt] <- shortest.paths(graph=graph.adjacency(Y[,,tt, 1]),mode="all")
  for(tt in 1:TT) dissim[,,tt] <- igraph::shortest.paths(graph=igraph::graph.adjacency(apply(Y[ , , tt, ], c(1, 2), mean)),mode="all") # sum over K (doesn't help), just doing by one time point for now
  dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
  X <- list()
  X[[1]] <- array(0,c(d,TT,n))
  X[[1]][,1,] <- t(cmdscale(d=dissim[,,1],k=d))
  temp.lambda <- 10
  H <- matrix(-1/n,n,n)+diag(n)

  for(tt in 2:TT){
    temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
      temp.lambda/(1+temp.lambda)*t(X[[1]][,tt-1,])%*%X[[1]][,tt-1,]
    temp <- eigen(temp)
    X[[1]][,tt,] <- t(temp$vectors[,1:d]%*%diag(temp$values[1:d])^(1/2))
    X[[1]][,tt,] <- t(vegan::procrustes(X=t(X[[1]][,tt-1,]),Y=t(X[[1]][,tt,]),scale=FALSE)$Yrot)
  }

  dim(X[[1]])

  initialize.wrap <- function(x){
    -c_initialize1(X[[1]],
                   c(n,d,TT),
                   Y,XSCALE=1/n,
                   BETAIN=x[1],BETAOUT=x[2],w[,1])
  }
  initialize.grad.wrap <- function(x){
    -c_initialize1_grad(X[[1]],
                        c(n,d,TT),
                        Y,XSCALE=1/n,
                        BETAIN=x[1],BETAOUT=x[2],w[,1])[2:3]
  }
  Optim <- optim(par=c(1,1),fn=initialize.wrap,
                 gr=initialize.grad.wrap,method="BFGS")


  X[[1]] <- X[[1]]/n
  # Bin[1] <- max(Optim$par[1],1e-4)
  # Bout[1] <- max(Optim$par[2],1e-4)

  Xit0 <- t(X[[1]][,1,])
  for(tt in 2:TT)Xit0 <- rbind(Xit0, t(X[[1]][,tt,]))
  Xit0 <- Xit0 -
    kronecker(rep(1,n*TT),matrix(apply(Xit0,2,mean),1,d))

  Xit0 <- scale(Xit0)

  # Brian:
  t <- 1
  Z_0_sewell <- array(NA, dim = c(n, 2, TT))
  for(t in 1:tt) {
    index_t <- (t-1)*n # batch of n nodes at time t
    Z_0_sewell[ , , t] <- Xit0[(index_t + 1):(index_t + n), ]
  }
  Z_0_sewell <- Z_0_sewell * sqrt(sigma2_0)

  return(Z_0_sewell)
}


#' Used for initialising alpha, see Gwee's paper and code where he uses this idea
#'
#' @param y
#' @param x
#' @param intercept
#' @param mode
#' @param diag
#'
#' @return
#' @export
#'
#' @examples
netpois <- function (y, x, intercept = TRUE, mode = "digraph", diag = FALSE)
{
  gfit <- function(glist, mode, diag) {
    y <- sna::gvectorize(glist[[1]], mode = mode, diag = diag,
                    censor.as.na = TRUE)
    x <- vector()
    for (i in 2:length(glist)) x <- cbind(x, sna::gvectorize(glist[[i]],
                                                        mode = mode, diag = diag, censor.as.na = TRUE))
    if (!is.matrix(x))
      x <- matrix(x, ncol = 1)
    mis <- is.na(y) | apply(is.na(x), 1, any)
    glm.fit(x[!mis, ], y[!mis], family = poisson(), intercept = FALSE)
  }
  if (is.list(y) || ((length(dim(y)) > 2) && (dim(y)[1] >
                                              1)))
    stop("y must be a single graph in netlogit.")
  if (length(dim(y)) > 2)
    y <- y[1, , ]
  if (is.list(x) || (dim(x)[2] != dim(y)[2]))
    stop("Homogeneous graph orders required in netlogit.")
  nx <- sna::stackcount(x) + intercept
  n <- dim(y)[2]
  g <- list(y)
  if (intercept)
    g[[2]] <- matrix(1, n, n)
  if (nx - intercept == 1)
    g[[2 + intercept]] <- x
  else for (i in 1:(nx - intercept)) g[[i + 1 + intercept]] <- x[i,
                                                                 , ]
  if (any(sapply(lapply(g, is.na), any)))
    warning("Missing data supplied to netlogit; this may pose problems for certain null hypotheses.  Hope you know what you're doing....")
  fit.base <- gfit(g, mode = mode, diag = diag)
  fit <- list()
  fit.base$coefficients
}

#' Initialise the alphas using the regression method
#'
#' @param Y
#' @param Z
#'
#' @return
#' @export
#'
#' @examples
initialise_alpha <- function(Y, Z) {
  alpha_0 <- matrix(NA, nrow = K, ncol = TT)
  for(t in 1:TT) {
    for(k in 1:K) {
      alpha_0[k, t] <- netpois(Y[ , , t, k], Dist(Z[ , , t]))[1]
    }
  }
  alpha_0
}
