devtools::load_all(path = 'other-code/simsurv')
library(tibble)
library(dplyr)
library(magrittr)
library(GGally)
library(sn)

set.seed(12931)

generate_prior_sample <- function() {
  res <- list(
    hazard_gamma = rgamma(n = 1, shape = 9.05, rate = 8.72),
    w_i = c(1, rnorm(n = 18)),
    theta = sn::rsn(n = 19, xi = 0, omega = 0.5, alpha = -1),
    alpha = sn::rsn(n = 1, xi = 0, omega = 0.5, alpha = -2),
    k_i = abs(rnorm(n = 1, mean = 4.5)),
    eta_before = abs(rnorm(n = 1, mean = 5, sd = 2.5)),
    eta_after = abs(rnorm(n = 1, mean = 5, sd = 2.5)),
    dd_gamma = rgamma(n = 1, shape = 2, rate = 2)
  )
  
  res$baseline_term <- (res$w_i %*% res$theta) %>% as.numeric()
  return(res)
}

pointwise_unnormalised_log_hazard <- function(t, x, betas) {
  with(betas, {
    log_hazard_term <- log(hazard_gamma) + ((hazard_gamma - 1) * log(t)) + (baseline_term)
    
    if (t < k_i) {
      longitudinal_term <- alpha * eta_before
      log_hazard_term <- log_hazard_term + longitudinal_term
    } else {
      longitudianl_term_2 <- alpha * eta_after
      log_hazard_term <- log_hazard_term + longitudianl_term_2
    }
    
    log_dd_term <- log(dd_gamma)
    log_res <- matrixStats::logSumExp(c(log_hazard_term, log_dd_term))
    
    return(log_res)
  })
}

n_sims <- 20000

many_df <- lapply(1 : n_sims, function(x) {
  generate_prior_sample() %>% 
    inset2('w_i', NULL) %>% 
    inset2('theta', NULL) %>% 
    as_tibble()
})

res <- parallel::mclapply(many_df, mc.cores = 6, function(z) {
  tryCatch(
    expr = {
      simsurv(
        loghazard = pointwise_unnormalised_log_hazard,
        x = tibble(id = runif(1)),
        betas = z,
        interval = c(1e-12, 200),
        rootsolver = 'uniroot',
        seed = 12931
      )
    },
    error = function(e) {
      e
    }
  )
})

samples_with_errors <- lapply(1 : n_sims, function(sample_id) {
  if ('error' %in% class(res[[sample_id]])) {
    event_time <- NA
    simsurv_status <- res[[sample_id]][['message']]
  } else {
    event_time <- res[[sample_id]][['eventtime']]
    simsurv_status <- 'Ok'
  }
  
  inner_res <- many_df[[sample_id]]
  inner_res$event_time <- event_time
  inner_res$simsurv_status <- simsurv_status
  
  return(inner_res)
  
}) %>% 
  bind_rows() %>%
  mutate(id = 1 : n())

ggpairs(
  data = samples_with_errors %>% 
    filter(is.na(event_time) | id < 5000),
  columns = 1 : 8,
  mapping = aes(colour = simsurv_status),
  upper = list(continuous = function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_hex(bins = 20)
  }),
  lower = list(continuous = function(data, mapping, ...) {
    ggally_autopoint(data, mapping, alpha = 0.1)
  })
)

plot_t_vals <- seq(from = -15, to = 5, length.out = 500) %>% 
  exp()

plot_ids <- samples_with_errors %>% 
  filter(simsurv_status == 'Ok') %>% 
  slice_sample(n = 14) %>% 
  pull(id) %>% 
  c(which(is.na(samples_with_errors$event_time)))

plot_tbl <- parallel::mclapply(plot_ids, mc.cores = 6, function(iter_id) {
  log_hazard_values <- array(NA, dim = length(plot_t_vals))
  input_df <- samples_with_errors[iter_id, ] %>% 
    as.data.frame()
  
  for (ii in seq_len(length(plot_t_vals))){
    log_hazard_values[ii] <- pointwise_unnormalised_log_hazard(
      t = plot_t_vals[ii],
      x = NULL,
      betas = as.list(input_df)
    )
  }
  
  tibble(
    id = iter_id,
    plot_t = plot_t_vals,
    log_hazard_values = log_hazard_values,
    label_str = with(input_df, 
      sprintf(
        'haz_gamma=%.2f, alpha=%.2f, eta_b=%.2f, eta_a=%.2f, base=%.2f',
        hazard_gamma,
        alpha,
        eta_before,
        eta_after,
        baseline_term
      )
    )
  )
}) %>% 
  bind_rows() %>% 
  left_join(samples_with_errors, by = 'id')

ggplot(
  plot_tbl,
  aes(x = plot_t, y = log_hazard_values, col = simsurv_status)
) +
  geom_line() +
  facet_wrap(vars(label_str))

