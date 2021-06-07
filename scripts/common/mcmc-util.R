library(coda)
library(bayesplot)
library(ggplot2)
library(knitr)
library(kableExtra)
library(stringr)
library(patchwork)
library(dplyr)

mcmc_list_to_array <- function(mcmc_list) {
  stopifnot(class(mcmc_list) == "mcmc.list")
  
  n_iter <- nrow(mcmc_list[[1]])
  n_par <- ncol(mcmc_list[[1]])
  n_chain <- length(mcmc_list)
  par_names <- dimnames(mcmc_list[[1]])[[2]]
  
  results_array <- array(
    data = NA,
    dim = c(n_iter, n_chain, n_par),
    dimnames = list(
      NULL,
      sprintf("chain_%d", 1:n_chain),
      par_names
    )
  )
  
  for (chain_index in 1:n_chain) {
    results_array[, chain_index,] <- mcmc_list[[chain_index]]
  }
  
  return(results_array)
  
}

array_to_mcmc_list <- function(array) {
  stopifnot(length(dim(array)) == 3)

  n_chain <- dim(array)[2]
  res <- lapply(1 : n_chain, function(chain_id) {
    coda::mcmc(array[ , chain_id, ])
  })

  return(coda::as.mcmc.list(res))
}


# takes an Iter X chain X param array returns a ggplot
plot_worst_pars <- function(
  samples_array,
  facet_name_value_pairs, # named vector pair, e.g. c('event_time' = 'italic(t)')
  n_warmup = 1
) {
  numerical_diags <- rstan::monitor(samples_array, n_warmup, print = FALSE) %>%
    as.data.frame()

  worst_index <- c(
    which.max(numerical_diags$Rhat),
    which.min(numerical_diags$n_eff)
  ) %>%
    unique()
  
  if (length(worst_index) == 1) {
    # if the worst rhat and neff are the same variable, 
    # get the second worst rhat
    worst_index <- c(
      worst_index,
      order(numerical_diags$Rhat, decreasing = TRUE)[2]  
    )
  }
  
  current_names <- names(samples_array[1, 1, worst_index])
  n_samples <- dim(samples_array)[1]
  
  # read str_replace_all doc VERY CAREFULLY! Caught out by this behaviour
  # for the second time now.
  # the named vector behaviour of str_replace_all all is inconsistent
  # I suggest passing a closure through the function argument (more consistent)
  # and matching on ".+" -- any char vec, which is more consistent
  if (is.function(facet_name_value_pairs)) {
    ideal_names <- str_replace_all(
      current_names,
      ".+",
      replacement = facet_name_value_pairs
    )
  } else {
    ideal_names <- str_replace_all(current_names, facet_name_value_pairs)    
  }
  
  names(ideal_names) <- current_names
  my_lab <- as_labeller(ideal_names, label_parsed)

  trace <- mcmc_trace(samples_array[(n_warmup + 1) : n_samples, , worst_index]) + 
    facet_wrap("parameter", ncol = 1, scales = "free_y", labeller = my_lab) + 
    xlab("Iteration") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = rel(0.8)),
      axis.text.y = element_text(size = rel(0.8))
    ) +
    bayesplot:::force_x_axis_in_facets()
  
  rank <- mcmc_rank_overlay(samples_array[, , worst_index], ref_line = TRUE) + 
    facet_wrap("parameter", ncol = 1, labeller = my_lab) + 
    theme(
      axis.text.x = element_text(size = rel(0.8)),
      axis.text.y = element_text(size = rel(0.8))
    )
  
  trace + rank + plot_layout(guides = "collect")
}

# takes the same array as above, and some extra things, and writes
# the appropriate latex table to a .tex file.
write_diagnostics_to_file <- function(
  samples_array,
  row_latex_names,
  output_filename,
  monitor_cols = c("mean", "Q5", "Q95", "Rhat", "Bulk_ESS", "Tail_ESS"),
  col_latex_names = c(
    "Mean",
    r"($q_{0.05}$)",
    r"($q_{0.95}$)", 
    r"($\widehat{R}$)", 
    r"($\widehat{\text{ESS}}_{B}$)", 
    r"($\widehat{\text{ESS}}_{T}$)"
  ),
  n_warmup = 1,
  n_digits = 2,
  ...
) {
  table <- rstan::monitor(samples_array, n_warmup, print = FALSE) %>%
    as.data.frame() %>%
    select(!!!monitor_cols)
  
  rownames(table) <- row_latex_names
  formatted_kable <- kable(
    table, 
    format = "latex", 
    digits = n_digits, 
    col.names = col_latex_names,
    escape = FALSE,
    booktabs = T,
    longtable = T,
    linesep = "",
    ...
  ) %>%
    kable_styling(latex_options = c("striped", "repeat_header"))
  
  invisible(cat(formatted_kable, file = output_filename))
}

bind_named_sublists <- function(outer_list, name, along_dim) {
  lapply(outer_list, function(sub_list) sub_list[[name]]) %>%
    abind(along = along_dim)
}