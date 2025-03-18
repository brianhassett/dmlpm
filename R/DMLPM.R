#' Implements the DMLPM model
#'
#' @param num_samples
#' @param burn
#' @param thin
#' @param post_pred_thin_factor
#' @param num_chains
#' @param Y
#' @param sigma2_0
#' @param sigma2_a
#' @param sigma2_z
#' @param alpha_0
#' @param Z_0
#' @param eta_alpha
#' @param eta_z
#' @param tune
#' @param delta_alpha
#' @param delta_z
#' @param AR_min
#' @param AR_max
#' @param discard
#'
#' @return
#' @export
#'
#' @examples
DMLPM <- function(num_samples = 2e3, burn = 10e3,
                  thin = 10,                        # store posterior sample every `thin` iterations
                  post_pred_thin_factor = 5,        # sample post. pred. `num_samples/post_pred_thin_factor` times
                  num_chains = 1,                   # no. of MCMC chains (parallel via the foreach package)
                  Y,                                # data
                  sigma2_0 = sigma2_0, sigma2_a = 1, sigma2_z = 1, # priors
                  alpha_0 = 1, Z_0 = NULL,          # initial values
                  eta_alpha = 1, eta_z = 1,         # std. dev. of random walk proposals
                  tune = TRUE,                      # turn on tuning of proposal SD's
                  delta_alpha = 0.1, delta_z = 0.1, # delta is the proportional change in SD during tuning
                  AR_min = 0.2, AR_max = 0.4,       # AR_min and AR_max provide a target range for the AR
                  discard = FALSE) {                # set to TRUE to discard burn-in and thinned out samples

  # Start timer:
  start1 <- start_timer()

  # No. of nodes, layers and time points:
  n <- nrow(Y)
  d <- 2 # change this later to an argument if I want to exlplore different niumbers of dimensions
  if(is.matrix(Y) | dim(Y)[3] == 1) { # single view
    K <- 1  # K = 1
    TT <- 1 # T = 1
  } else { # multiview
    K <- dim(Y)[4]  # K layers
    TT <- dim(Y)[3] # T time points
  }

  # No. of iterations:
  # (total no. of iterations required for `num_samples` posterior samples with thinning)
  iters <- num_samples*thin

  # No. of posterior predictive samples:
  # (combine thin and post_pred_thin_factor to compute this)
  num_post_pred_samples <- num_samples/post_pred_thin_factor
  post_pred_thin <- thin*post_pred_thin_factor

  # Parallel:
  if(num_chains > 1) {
    library(doMC)
    registerDoMC(cores = num_chains)
  }

  mcmc_chains <- foreach(chain_id = 1:num_chains) %dopar% {

    # Create log file:
    log_file <- paste0("logs/chain", chain_id, ".txt")
    file.create(log_file)

    # Arrays, matrices, and vectors to store output:
    log_posterior <- vector("double", length = (burn + iters))
    log_posterior_thin = vector("double", length = num_samples)
    Z_out <- array(NA, dim = c(n, d, TT, (burn + iters)))
    Z_out_thin = array(NA, dim = c(n, d, TT, num_samples))
    Z_acc <- array(NA, dim = c(n, TT, (burn + iters))) # no need for the 2 columns for each dim, since both dims are accepted/rejected together anyway, so this is a 3D array like with the multiview version
    alpha_out <- array(NA, dim = c(K, TT, (burn + iters))) # alpha is now a K x T matrix, so stack one of these for each iter along 3rd dim
    alpha_out_thin = array(NA, dim = c(K, TT, num_samples))
    alpha_acc <- array(NA, dim = c(K, TT, (burn + iters)))
    # alpha_out <- matrix(NA, nrow = (burn + iters), ncol = K)
    # alpha_out_thin = matrix(NA, nrow = (num_samples), ncol = K)
    # alpha_acc <- matrix(NA, nrow = (burn + iters), ncol = K)

    # Initialise latent positions if not provided:
    # if(is.null(Z_0))
    #   Z_0 <- mds_initialise(Y[ , , 1])

    # For multiple chains add perturbation to initial values:
    # if(num_chains > 1) {
    #   epsilon <- min(abs(Z_0))/10 # perturbation SD (keep this small to avoid huge relative perturbations for some nodes)
    #   Z_0 <- Z_0 + rnorm(n*ncol(Z_0), mean = 0, sd = epsilon)
    #   alpha_0 <- alpha_0 + rnorm(TT, mean = 1, sd = 0.1)
    # }

    # Store initial values:
    log_posterior[1] <- log_posterior(Y = Y, Z = Z_0, alpha = alpha_0, sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z)
    log_posterior_cur <- log_posterior[1]
    Z_out[, , , 1] <- Z_0  # -same initial positions for all t     # use different initial values for each chain
    Z_acc[, , 1] <- 0         # reject initial value (part of burn-in anyway)
    alpha_out[ , , 1] <- alpha_0 # use different initial values for each chain
    alpha_acc[ , , 1] <- 0       # reject initial value (part of burn-in anyway)

    # Posterior predictive checks:
    # (for the dynamic model the t subscript here means iter, not time, this is confusing so I will change it)
    # For the dynamic multiview model I will change p_ijt to just P_hat because it will be too confusing with all the subscripts.
    # Another (4th) dimension needs to be added for the layers since the p_ij depends on time AND layer (alpha_k)
    #L_hat <- array(NA, dim = c(num_post_pred_samples, choose(n, 2), TT, K)) # each row is the pairwise prob's based on the current (iter's) Z values, the prob's corresponding to each time point are stacked along the 3rd dimension, the prob's corresponding to each layer k are stacked along the 4th dimension
    #Y_hat <- array(NA, dim = c(num_post_pred_samples, choose(n, 2), TT, K)) # each row is sampled from Bernoulli distributions with prob's given by P_hat, the samples for each time point are stacked along the 3rd dimension, the samples for each layer k are stacked along the 4th dimension

    # Proposal SDs:
    # (these are objects to store the proposal SDs as they are updated,
    #  the acceptance rate is checked after every 1000 burn-in iterations)
    eta_alpha_i <- array(NA, dim = c(K, TT, burn/1000 + 1)) # the proposal SDs for the alphas, K rows, T cols, with the updated SD's along 3rd dim
    eta_alpha_i[ , , 1] <- eta_alpha # store initial
    # eta_alpha_i <- matrix(NA, nrow = K, ncol = burn/1000 + 1) # the proposal SDs for the alphas, K rows with the updated SD's along the cols
    # eta_alpha_i[ , 1] <- eta_alpha # store initial
    AR_alpha_i <- array(NA, dim = c(K, TT, burn/1000)) # AR for each batch of 1000 samples (one less entry since we compute AR after first 1000)
    eta_alpha <- eta_alpha_i[ , , 1] # the current proposal SD's are the K values in the first slice
    # AR_alpha_i <- matrix(NA, nrow = K, ncol = burn/1000) # AR for each batch of 1000 samples (one less entry since we compute AR after first 1000)
    # eta_alpha <- eta_alpha_i[ , 1] # the current proposal SD's are the K values in the first column
    eta_z_i <- array(NA, dim = c(n, TT, burn/1000 + 1)) # nodes on the rows, updated proposal SD on the columns
    eta_z_i[ , , 1] <- eta_z # store initial (same for all nodes and time points at the start)
    AR_z_i <- array(NA, dim = c(n, TT, burn/1000)) # AR for each batch of 1000 samples (one less column since we compute AR after first 1000)
    eta_z <- eta_z_i[ , , 1] # This will hold the current proposal SDs and will be updated

    write("Burning in ...", file = log_file, append = TRUE)
    # Sample and update using M-H:
    for (i in 2:(burn + iters)) { # TODO maybe change i to iter, frees up i for use later
      # Print progress
      if(i %% 10 == 0) {
        print(paste0("i = ", i, " (", Sys.time(), ")"))
      }
      if(i %% 1000 == 0) {
        line <- paste0("Chain ", chain_id, ": iteration ", i, "/", burn+iters)
        write(line, file = log_file, append = TRUE)
        if(i == burn) {
          write("... burn-in complete.", file = log_file, append = TRUE)
          write("Sampling ...", file = log_file, append = TRUE)
        }
      }

      # Current values:
      Z_cur <- Z_out[ , , , i - 1]
      alpha_cur <- alpha_out[ , , i - 1]

      # Oct 2024: (need this if burn == 0 to avoid NA's at the start of Z_out_thin):
      if(burn == 0)
        Z_out_thin[ , , , 1] <- Z_cur

      ###
      ### Update alpha ###
      ### (now need to update for each k and t)
      for(t in 1:TT) { # Do t first, so that I can compute D_t and reuse it for each t

        D_t <- Rfast::Dist(Z_cur[ , , t])

        for(k in 1:K) {

          alpha_kt_cur <- alpha_cur[k, t] # current values of alpha
          alpha_kt_prop <- rnorm(n = 1, mean = alpha_kt_cur, sd = eta_alpha[k, t]) # proposed value

          # Auxiliary value from U(0, 1):
          u <- runif(n = 1, min = 0, max = 1)

          # Compute log-posterior and log-acceptance function:
          alpha_prop <- alpha_cur
          alpha_prop[k, t] <- alpha_kt_prop
          #log_posterior_prop <- log_posterior(Y = Y, Z = Z_cur, alpha = alpha_prop, sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z)
          log_posterior_cur_kt <- log_posterior_alpha_kt(Y_tk = Y_tk_list[[t]][[k]],
                                                         D_t = D_t,
                                                         alpha = alpha_cur,
                                                         sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z,
                                                         k = k, t = t)
          log_posterior_prop_kt <- log_posterior_alpha_kt(Y_tk = Y_tk_list[[t]][[k]],
                                                          D_t = D_t,
                                                          alpha = alpha_prop,
                                                          sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z,
                                                          k = k, t = t)


          log_acc <- min(0, log_posterior_prop_kt - log_posterior_cur_kt) # for alpha[k, t]

          # Acceptance:
          if (log(u) < log_acc) {
            alpha_cur <- alpha_prop
            log_posterior_cur <- log_posterior_cur + log_posterior_prop_kt - log_posterior_cur_kt
            alpha_acc[k, t, i] <- 1 # accepted
          } else {
            alpha_acc[k, t, i] <- 0 # rejected (no need to update alpha_cur or log_posterior_cur)
          }
        }
      }

      # Store values: (Dec 7 2024: moved this outside the previous loop, as I think it makes more sense, otherwise I rewrite for each k and keep the last, makes sense to just store it once after the loop is finished)
      # It makes sense since all of the alpha_kt values have been updated, so now I can store the slice, otherwise I am just overwriting every time.
      # This might actually save a bit of time as all this accessing could be slow,
      # so change it in the other code when I am combining everything.
      alpha_out[ , , i] <- alpha_cur
      if(i > burn & i %% thin == 0) { # store thinned values
        alpha_out_thin[ , , (i - burn)/thin] <- alpha_cur
      }

      ###
      ### Update Z ###
      ###
      for(t in 1:TT) {

        for (j in 1:n) {
          # Current z_j value:
          z_j_cur <- c(Z_cur[j, 1, t], Z_cur[j, 2, t]) # (z_j1, z_j2), i.e. the position of the jth node at time t

          # Proposed z_j value:
          # (just use two independent draws from univariate normals)
          z_j1_prop <- rnorm(n = 1, mean = z_j_cur[1], sd = eta_z[j, t])
          z_j2_prop <- rnorm(n = 1, mean = z_j_cur[2], sd = eta_z[j, t])
          z_j_prop <- c(z_j1_prop, z_j2_prop)

          # Auxiliary value from U(0, 1):
          u <- runif(1, 0, 1)

          # Compute log-posterior and log-acceptance function:
          Z_prop <- Z_cur
          Z_prop[j, , t] <- z_j_prop
          log_posterior_cur_j <- log_posterior_it_cpp(Y = Y_list, Z = Z_cur, alpha = alpha_cur, i = j - 1, t = t - 1, sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z) # subtract 1 for zero-indexing

          log_posterior_prop_j <- log_posterior_it_cpp(Y = Y_list, Z = Z_prop, alpha = alpha_cur, i = j - 1, t = t - 1, sigma2_0 = sigma2_0, sigma2_a = sigma2_a, sigma2_z = sigma2_z) # subtract 1 for zero-indexing

          log_acc <- min(0, log_posterior_prop_j - log_posterior_cur_j) # for the ith node
          #log_acc <- min(0, log_posterior_prop - log_posterior_cur)

          # Acceptance:
          if (log(u) < log_acc) {
            z_j_cur <- z_j_prop
            log_posterior_cur <- log_posterior_cur + log_posterior_prop_j - log_posterior_cur_j # update the full posterior using the node i info
            Z_acc[j, t, i] <- 1 # accepted
          } else {
            Z_acc[j, t, i] <- 0 # rejected
          }

          # Update values in Z_cur before moving on to the next node:
          Z_cur[j, , t] <- z_j_cur
        }
      }


      # Store the i-th Z matrix (now that all rows (nodes) have been updated):
      Z_out[, , , i] <- Z_cur
      if(i > burn & i %% thin == 0) {
        Z_out_thin[ , , , (i - burn)/thin] <- Z_cur # store thinned values after a full update of Z
      }

      # Store the current log-posterior after each iteration i:
      log_posterior[i] <- log_posterior_cur
      if(i > burn & (i - burn) %% thin == 0) # store thinned values
        log_posterior_thin[(i - burn)/thin] <- log_posterior_cur

      # Oct 2024: changed this for the Poisson mode - just needed to compute lambdas instead of prob's, and then simulate random Poisson data.
      # Posterior predictive:
      # (for dynamic multiview I just added an extra (outer) loop for the K layers and indexed in to the arrays and alpha vector appropriately)
      #if(i > burn & (i - burn) %% post_pred_thin == 0) { # don't think I need the `i - burn` part, should just be the same as normal thinning, see above
      #  for(k in 1:K) {
      #    for(t in 1:TT) {
      #      # P_hat[(i - burn)/post_pred_thin, , t, k] <- pairwise_probabilities_vector(alpha = alpha_cur[k], Z = Z_cur[ , , t]) # vector of pairwise probabilities
      #      # Y_hat[(i - burn)/post_pred_thin, , t, k] <- sapply(P_hat[(i - burn)/post_pred_thin, , t, k], function(p) sample(0:1, size = 1, prob = c(1-p, p)))
      #      L_hat[(i - burn)/post_pred_thin, , t, k] <- pairwise_lambda_vector(alpha = alpha_cur[k], Z = Z_cur[ , , t]) # vector of pairwise lambdas (rates)
      #      Y_hat[(i - burn)/post_pred_thin, , t, k] <- sapply(L_hat[(i - burn)/post_pred_thin, , t, k], function(lambda) rpois(n = 1, lambda))
      #    }
      #   }
      # }

      ### Tune variance (SD) every 1,000 iterations during burn-in ###
      if(tune) {
        if ( (i %% 1000 == 0) & (burn >= i) ) {

          # Compute and store AR for the last 1000 iterations:
          # alpha
          AR_alpha <- apply(alpha_acc[ , , (i - 999):i], MARGIN = c(1, 2), mean)
          AR_alpha_i[ , , i/1000] <- AR_alpha
          # AR_alpha <- colMeans(as.matrix(alpha_acc[(i - 999):i, ])) # convert back to matrix when K = 1
          # AR_alpha_i[ , i/1000] <- AR_alpha

          # Z
          AR_Z <- matrix(NA, nrow = n, ncol = TT) # nodes on rows, time points on columns
          for(t in 1:TT) {
            for(k in 1:n) {
              AR_Z[k, t] <- mean(Z_acc[k, t, (i - 999):i])
            }
          }
          AR_z_i[ , , i/1000] <- AR_Z

          # Tune proposal variance (SD) for alpha:
          for(k in 1:K) {
            for(t in 1:TT) {
              if( AR_alpha[k, t] < AR_min ) {
                # AR is too low, decrease variance:
                eta_alpha[k, t] <- eta_alpha[k, t] - delta_alpha*eta_alpha[k, t] # update
                eta_alpha_i[k, t, i/1000 + 1] <- eta_alpha[k, t]              # store
              } else if( AR_alpha[k, t] > AR_max ) {
                # AR is too high, increase variance:
                eta_alpha[k, t] <- eta_alpha[k, t] + delta_alpha*eta_alpha[k, t] # update
                eta_alpha_i[k, t, i/1000 + 1] <- eta_alpha[k, t]              # store
              } else {
                # AR is within limits, maintain the variance and store it:
                eta_alpha_i[k, t, i/1000 + 1] <- eta_alpha[k, t] # store
              }
            }
          }

          # Tune proposal variance for Z:
          for(t in 1:TT) {
            for(k in 1:n) {
              if( AR_Z[k, t] < AR_min ) {
                # AR is too low, decrease variance:
                eta_z_i[k, t, i/1000 + 1]  <- eta_z_i[k, t, i/1000] - delta_z*eta_z_i[k, t, i/1000] # update and store
              } else if( AR_Z[k, t] > AR_max ) {
                # AR is too high, increase variance:
                eta_z_i[k, t, i/1000 + 1]  <- eta_z_i[k, t, i/1000] + delta_z*eta_z_i[k, t, i/1000] # update and store
              } else {
                # AR is within limits, maintain the variance and store it:
                eta_z_i[k, t, i/1000 + 1] <- eta_z_i[k, t, i/1000] # store
              }
            }
          }

          # Update the current vector of proposal SDs for Z:
          # (the alpha SD's are already updated on the fly in eta_alpha)
          eta_z <- eta_z_i[ , , i/1000 + 1]
        }
      }
    }

    # Compute mean acceptance rates:
    alpha_acc_mean <- apply(alpha_acc[, , (burn + 1):(iters + burn)], c(1, 2), mean)
    Z_acc_mean <- apply(Z_acc[, , (burn + 1):(iters + burn)], c(1, 2), mean)

    # Old loops:
    #alpha_acc_mean <- matrix(NA, nrow = K, ncol = TT) #<- colMeans(alpha_acc)
    # for(i in 1:n) {
    #   for(t in 1:TT) {
    #     alpha_acc_mean[i, t] <- mean(Z_acc[i, t, (burn + 1):(iters + burn)])
    #   }
    # }

    #Z_acc_mean <- matrix(NA, nrow = n, ncol = TT)

    # for(i in 1:n) {
    #   for(t in 1:TT) {
    #     Z_acc_mean[i, t] <- mean(Z_acc[i, t, (burn + 1):(iters + burn)])
    #   }
    # }

    write("... sampling complete.", file = log_file, append = TRUE)

    # Return the output as a list:
    if (discard) { # discard burn-in and thinned out samples
      return(list(
        iters = iters, burn = burn,
        thin = thin, post_pred_thin = post_pred_thin,
        num_samples = num_samples, num_post_pred_samples = num_post_pred_samples,
        num_chains = num_chains, chain_id = chain_id,
        n = n, TT = TT, K = K,
        log_posterior_thin = log_posterior_thin,
        alpha_out_thin = alpha_out_thin,
        alpha_acc_mean = alpha_acc_mean,
        Z_out_thin = Z_out_thin,
        Z_acc_mean = Z_acc_mean,
        eta_alpha_i = eta_alpha_i, eta_z_i = eta_z_i,
        delta_alpha = delta_alpha, delta_z = delta_z,
        AR_min = AR_min, AR_max = AR_max#,
        #L_hat = L_hat,
        #Y_hat = Y_hat
      ))
    } else { # keep burn-in and all posterior samples
      return(list(
        iters = iters, burn = burn,
        thin = thin, post_pred_thin = post_pred_thin,
        num_samples = num_samples, num_post_pred_samples = num_post_pred_samples,
        num_chains = num_chains, chain_id = chain_id,
        n = n, TT = TT, K = K,
        log_posterior = log_posterior,
        log_posterior_thin = log_posterior_thin,
        alpha_out = alpha_out,
        alpha_acc = alpha_acc, alpha_acc_mean = alpha_acc_mean,
        alpha_out_thin = alpha_out_thin,
        Z_out = Z_out,
        Z_out_thin = Z_out_thin,
        Z_acc = Z_acc, Z_acc_mean = Z_acc_mean,
        eta_alpha_i = eta_alpha_i, eta_z_i = eta_z_i,
        delta_alpha = delta_alpha, delta_z = delta_z,
        AR_min = AR_min, AR_max = AR_max#,
       # L_hat = L_hat,
        #Y_hat = Y_hat
      ))
    }
  }

  ###
  ### Post processing
  ### (for now this just uses the MAP chain to compute the posterior means.
  ###  might update this later.)
  ###

  # MAP estimate:
  i_map <- data.frame(chain = 1:num_chains,
                      i = map_int(mcmc_chains, function(chain) which.max(chain$log_posterior_thin)),
                      log_posterior_max = map_dbl(mcmc_chains, function(chain) max(chain$log_posterior_thin)))

  # Store the MAP estimates of alpha and Z:
  lp_max <- slice_max(i_map, order_by = log_posterior_max)
  alpha_map <- mcmc_chains[[lp_max$chain]]$alpha_out_thin[ , , lp_max$i]
  Z_map <- mcmc_chains[[lp_max$chain]]$Z_out_thin[ , , , lp_max$i]

  # Just one chain for now:
  chain <- lp_max$chain

  ### Procrustes transformation with stacked nT x 2 matrices ###
  ###                                                        ###

  # First, construct the nT x 2 x num_samples array of all the stacked matrices:
  Z_stack <- array(NA, dim = c(n*TT, 2, mcmc_chains[[chain]]$num_samples))
  for(t in 1:TT) {
    Z_stack[((t - 1)*n + 1):(t*n), , ] <- mcmc_chains[[chain]]$Z_out_thin[ , , t, ] # These indexes look weird but they are correct:
  }

  # Now extract the MAP estimate (same MAP code as above works here):
  Z_ref <- Z_stack[ , , lp_max$i] # this is the nT x 2 reference matrix

  # Transform:
  Z_trans <- array(NA, dim = dim(Z_stack)) # Array to store the transformed samples dim = (nT, 2, num_samples)
  for(i in 1:dim(Z_trans)[3]) {
    Z_trans[ , , i] <- vegan::procrustes(X = Z_ref, Y = Z_stack[ , , i], scale = FALSE)$Yrot
  }

  # Unstack the matrices again:
  Z_trans_t <- array(NA, dim = dim(mcmc_chains[[chain]]$Z_out_thin))
  for(t in 1:TT) {
    Z_trans_t[ , , t, ] <- Z_trans[((t - 1)*n + 1):(t*n), , ]
  }

  # Posterior means for each Zt:
  Z_hat_t <- array(NA, dim = dim(Z_trans_t[ , , , 1]))
  for(t in 1:TT) {
    for(i in 1:n) {
      for(j in 1:2) {
        Z_hat_t[i, j, t] <- mean(Z_trans_t[i, j, t, ])
      }
    }
  }

  # Posterior mean for alpha:
  # Just the MAP chain:
  alpha_hat <- apply(mcmc_chains[[chain]]$alpha_out_thin, c(1, 2), mean)
  #alpha_hat <- colMeans(mcmc_chains[[chain]]$alpha_out_thin)

  # Stop timer and print run time:
  stop1 <- stop_timer()
  print(run_time <- print_run_time(start1, stop1))

  # Return:
  return(list(# MCMC:
    mcmc_chains = mcmc_chains,
    num_chains = num_chains,
    num_samples = num_samples,
    num_post_pred_samples = num_chains*iters/post_pred_thin,
    iters = iters, burn = burn,
    thin = thin, post_pred_thin_factor = post_pred_thin_factor,
    # MAP:
    map = list(chain_map = chain, i_map = i_map, lp_max = lp_max, alpha_map = alpha_map, Z_map = Z_map),
    # Data:
    Y = Y,                 # Data matrix/array
    n = n,                 # No. of nodes
    TT = TT,               # No. of time point
    K = K,                 # No. of layers
    # Priors:
    sigma2_a = sigma2_a,   # Prior variance on alpha
    sigma2_z = sigma2_z,   # Prior variance on z_ij
    # Posterior:
    Z_hat_t = Z_hat_t,     # Estimated posterior latent positions
    Z_trans = Z_trans,     # Stacked and transformed matrices for comparison later
    Z_trans_t = Z_trans_t,  # UnStacked and transformed matrices for comparison later
    alpha_hat = alpha_hat, # Estimated posterior intercept(s)
    # Misc:
    time = list(start = start1, stop = stop1, runtime = run_time),
    discard = discard))

  # Post-processing code goes here
  # Stop timer and print run time:
  stop1 <- stop_timer()
  print(run_time <- print_run_time(start1, stop1))
}
