#' Extract a ggnetwork data frame for an adjacency matrix
#'
#' @param adjacency_matrix
#' @param layout
#' @param include_unconnected
#'
#' @return
#' @export
#'
#' @examples
ggnetwork_df <- function(adjacency_matrix, layout, include_unconnected = FALSE) {
  # Set the node sizes (by Page-Rank) and colours:
  g1 <- network::as.network(adjacency_matrix)
  d <- degree(graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected"))
  g1 %v% "degree" <- d
  g1 %v% "color" <- ifelse(d > 0, "steelblue", "orange3")

  # Create a data frame from the network and its attributes:
  G1 <- ggnetwork(g1, scale = F, layout = layout)

  # This is a way to identify the nodes as opposed to edges:
  # (the ggnetwork function doesn't differentiate)
  G1 <- G1 %>%
    mutate(node_flag = near(x, xend) & near(y, yend))

  # Remove unconnected nodes from the plot:
  if(!include_unconnected) {
    G1 <- filter(G1, degree > 0)
  }

  return(G1)
}


#' Function to plot log-posterior:
#' (The default plot is a plot of the log-posterior of the thinned samples as I have been doing before.
#' I included some other options to plot all the samples, including burn-in or not,
#' and discarding a few values at the start to help with the skew of the scale.)
#'
#' @param fit
#' @param thinned
#' @param discard_burnin
#' @param exclude_initial
#' @param scales
#'
#' @return
#' @export
#'
#' @examples
plot_log_posterior <- function(fit, thinned = TRUE, discard_burnin = TRUE, exclude_initial = 0, scales = "fixed") {

  # Extract log-posterior for each chain:
  if (thinned) { # plot log-posterior of thinned samples
    df <- lapply(fit$mcmc_chains, function(chain) chain$log_posterior_thin) %>%
      set_names(paste0(1:fit$num_chains)) %>%
      bind_cols(i = 1:fit$num_samples) %>%
      pivot_longer(cols = !starts_with("i"), names_to = "chain", values_to = "log_posterior")
  } else { # plot log-posterior of all samples including burn-in
    df <- lapply(fit$mcmc_chains, function(chain) chain$log_posterior) %>%
      set_names(paste0(1:fit$num_chains)) %>%
      bind_cols(i = 1:(fit$burn + fit$iters)) %>%
      mutate(burn = c(rep(TRUE, fit$burn), rep(FALSE, fit$iters)))
    if(exclude_initial > 0) { # remove some initial values since they can skew the scale
      df <- slice_tail(df, n = -exclude_initial)
    }
    df <- df %>%
      pivot_longer(cols = (!starts_with("i") & !starts_with("burn")), names_to = "chain", values_to = "log_posterior")
    if(discard_burnin) {
      df <- filter(df, !burn)
    }
  }

  # Plot:
  P <- ggplot(df)

  if (thinned) { # plot log-posterior of thinned samples
    P <- P +
      geom_line(aes(i, log_posterior, colour = chain), linewidth = 0.25, alpha = 0.8, show.legend = F) +
      # geom_point(
      #   data = fit$map$i_map,
      #   aes(i, log_posterior_max), shape = "+", size = 3, colour = "black"
      # ) +
      # geom_point(
      #   data = slice_max(fit$map$i_map, log_posterior_max),
      #   aes(i, log_posterior_max), shape = 1, size = 4, colour = "black"
      # ) +
      facet_wrap(~chain, labeller = label_both, scales = scales) +
      #scale_color_viridis_d(direction = -1, option = "D", end = 0.95) +
      scale_color_manual(values = colour_scale) +
      labs(
        x = "i", y = "Log-posterior",
        title = paste0("Log-posterior after each (full) iteration (num_samples = ", fit$num_samples, ", thin = ", fit$thin, ")"),
        subtitle = "MAP for each chain marked with a cross, overall MAP marked with a circle"
      )
  } else { # plot log-posterior of all samples including burn-in
    P <- P +
      geom_line(aes(i, log_posterior, colour = burn), linewidth = 0.25, alpha = 0.8, show.legend = F) +
      scale_color_manual(values = c("steelblue", "grey")) +
      #bh_alpha_override_color() +
      facet_wrap(~chain, labeller = label_both, scales = scales) +
      labs(
        x = "i", y = "Log-posterior for all iterations (burn-in is in grey)",
        title = "Log-posterior",
      )
  }
  return(P)
}


## TODO
## Need to rework this function for the DM model, instead of just M or D.

#' Function to plot summary of posterior samples of alpha:
#' The facet argument applies to the trace plots only, the ACF must be faceted, the densities are fine since theya re easyo to see already
#'
#' @param fit
#' @param facet_trace
#'
#' @return
#' @export
#'
#' @examples
plot_alpha_posterior <- function(fit, facet_trace = FALSE) {

  # Posterior samples:
  alpha_df <- purrr::map(fit$mcmc_chains, function(chain) chain$alpha_out_thin) %>%
    map(data.frame) %>%
    bind_rows(.id = "chain") %>%
    mutate(i = rep(1:fit$num_samples, fit$num_chains), .before = 2) %>%
    mutate(chain = stringr::str_extract(chain, "\\d+")) %>%
    set_names(c("chain", "i", paste0("layer_", 1:fit$K))) %>%
    pivot_longer(cols = starts_with("layer"), names_to = "layer", values_to = "alpha") %>%
    mutate(layer = stringr::str_extract(layer, "\\d+") %>% as.numeric())
  alpha_true_df <- data.frame(layer = 1:fit$K, alpha = alpha)

  # Rhat for each layer:
  Rhat <- Rhat(fit$mcmc_chains, param = "alpha")

  # All ACF data (layers and chains):
  ACF <- map(split(alpha_df, alpha_df$layer), function(chains) {
    map(split(chains, chains$chain), function(chain) ggAcf(chain$alpha, plot = FALSE))
  }) %>%
    map(function(layer) map(layer, function(chain) data.frame(chain$lag, chain$acf) %>% set_names(c("lag", "acf"))) %>% bind_rows(.id = "chain")) %>%
    bind_rows(.id = "layer")

  # Plot trace, ACF, and densities for each layer (all chains):
  for (k in 1:fit$K) {
    # Trace plots:
    P1 <- ggplot() +
      geom_line(data = filter(alpha_df, layer == k), aes(i, alpha, group = chain, colour = chain), linewidth = 0.25, alpha = 0.7, show.legend = F) +
      geom_hline(data = filter(alpha_true_df, layer == k), aes(yintercept = alpha), colour = "black", linetype = 1) + # true
      geom_hline(yintercept = fit$alpha_hat[k], colour = "black", linetype = 2) + # posterior mean
      # labs(title = latex2exp::TeX(paste0("\\hat{$\\alpha$} = (", paste0(round(fit$alpha_hat, 2), collapse = ", ")))) # I had thought about adding the mean for each chain here but I haven't stored it, I accidentally sued the mean for each layer here, might come c=back to this later.
      labs(
        title = latex2exp::TeX(paste0("Trace plot for each chain ($\\hat{R} = ", format(round(Rhat[1, k], 2), nsmall = 2), "$)")),
        subtitle = latex2exp::TeX(paste0("$\\hat{\\alpha}_", k, "$ = ", round(fit$alpha_hat[k], 2), " (dashed line), $\\alpha_", k, "$ = ", round(alpha[k], 2), " (solid line)"))
      ) +
      scale_color_manual(values = colour_scale) +
      theme(plot.subtitle = element_text(size = 10))
    if(facet_trace)
      P1 <- P1 + facet_wrap(~chain)

    # ACF plots:
    ci <- qnorm((1 + 0.95) / 2) / sqrt(fit$num_samples)
    P2 <-
      ggplot(filter(ACF, layer == k, lag > 0)) +
      geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf, colour = chain), show.legend = F) +
      geom_hline(yintercept = c(-ci, ci), colour = "black", linetype = "solid") +
      facet_wrap(~chain, labeller = label_both) +
      labs(x = "Lag", y = "ACF", subtitle = paste0("ACF plot for each chain (thin = ", fit$thin, ")")) +
      scale_colour_manual(values = colour_scale) +
      theme(plot.subtitle = element_text(size = 12))

    # Density plots:
    P3 <- ggplot() +
      geom_density(data = filter(alpha_df, layer == k), aes(alpha, group = chain, colour = chain), alpha = 1, show.legend = F) +
      geom_vline(data = filter(alpha_true_df, layer == k), aes(xintercept = alpha), colour = "black", linetype = 1) + # true
      geom_vline(xintercept = fit$alpha_hat[k], colour = "black", linetype = 2) + # posterior mean
      # theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
      labs(
        title = "Posterior density for each chain",
        subtitle = latex2exp::TeX(paste0("$\\hat{\\alpha}_", k, "$ = ", round(fit$alpha_hat[k], 2), " (dashed line), $\\alpha_", k, "$ = ", round(alpha[k], 2), " (solid line)"))
      ) +
      scale_colour_manual(values = colour_scale) +
      theme(plot.subtitle = element_text(size = 10))
    print((P2 + P1 + P3) + plot_annotation(title = latex2exp::TeX(paste0("Layer ", k, ": $\\alpha_", k, "$ = ", alpha[k], "; $\\hat{\\alpha}_", k, "$ = ", round(fit$alpha_hat[k], 2))), subtitle = paste0("num_samples = ", fit$num_samples, "; burn = ", fit$burn, "; iters = ", fit$iters), theme = theme(plot.title = element_text(hjust = 0))))
  }
}
