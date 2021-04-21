library(dplyr)
library(tibble)
library(magrittr)

n_indiv <- 16

res <- lapply(1 : n_indiv, function(indiv_id) {
  obs_upper_limit <- rgamma(n = 1, shape = 5, scale = 3) %>% 
    round() %>% 
    subtract(0.5)    
  
  obs_lower_limit <- 0.5
  obs_times <- seq(
    from = obs_lower_limit, 
    to = obs_upper_limit,
    by = 1
  )
  
  n_obs <- length(obs_times)
  
  y_sigma <- rnorm(n = 1, mean = 50, sd = 25) %>% 
    abs()
  
  breakpoint <- runif(n = 1, min = obs_lower_limit, obs_upper_limit)
  intercept <- rnorm(n = 1, mean = 7500, sd = 1000) %>% 
    abs()
  
  slope_before <- rnorm(n = 1, mean = 2500, sd = 500) %>% 
    abs()
  
  slope_after <- rnorm(n = 1, mean = 1500, sd = 500) %>% 
    abs()
  
  obs_values <- array(dim = n_obs)
  
  for (ii in 1 : n_obs) {
    x_obs <- obs_times[ii]
    
    if (x_obs < breakpoint) {
      y_val <- intercept + slope_before * (x_obs - breakpoint)
    } else {
      y_val <- intercept + slope_after * (x_obs - breakpoint)
    }
      
    obs_values[ii] <- y_val
  }

  obs_values <- (obs_values + rnorm(n = n_obs, mean = 0, sd = y_sigma)) 
  
  sim_res <- tibble(
    indiv_id = indiv_id,
    n_obs = n_obs,
    obs_lower_limit = obs_lower_limit,
    obs_upper_limit = obs_upper_limit,
    y_sigma = y_sigma,
    breakpoint = breakpoint,
    intercept = intercept,
    slope_before = slope_before,
    slope_after = slope_after,
    x_obs = obs_times,
    y_obs = obs_values
  )
      
})

res <- res %>% 
  bind_rows() %>% 
  mutate(
    plot_label = sprintf(
      'i = %02d, beta_b = %.2f, beta_a = %.2f',
      indiv_id,
      slope_before,
      slope_after
    )
  )

library(ggplot2)

base_plot <- ggplot(res, aes(x = x_obs, y = y_obs)) +
  geom_point() +
  geom_vline(aes(xintercept = breakpoint)) +
  facet_wrap(vars(plot_label), scales = 'free') 

library(rstan)
library(purrr)

subset_vec <- res %>% 
  group_split(indiv_id) %>% 
  map_int(~ nrow(.)) %>% 
  cumsum() %>% 
  add(1) %>% 
  c(1, .)

breakpoint_limits <- res %>% 
  group_by(indiv_id) %>% 
  summarise(
    lower = min(x_obs),
    upper = max(x_obs)
  )

emp_prior_data <- res %>% 
  group_by(indiv_id) %>% 
  summarise(
    mean = mean(y_obs),
    sd = sd(y_obs)
  )

stan_data <- list(
  n_icu_stays = res %>% pull(indiv_id) %>% n_distinct(),
  n_total_obs = res %>% nrow(),
  subset_vector = subset_vec,
  y_vec = res %>% pull(y_obs), 
  x_vec = res %>% pull(x_obs),
  breakpoint_lower = breakpoint_limits %>% pull(lower),
  breakpoint_upper = breakpoint_limits %>% pull(upper),
  midpoint_prior_sd = 0.5
)

init_gen <- function(chain_id) {
  n_icu_stays <- stan_data$n_icu_stays
  breakpoint_midpoint <- (stan_data$breakpoint_lower + stan_data$breakpoint_upper) / 2
  
  list(
    breakpoint = array(
      data = rnorm(n = n_icu_stays, mean = breakpoint_midpoint, sd = 0.15) %>% abs(),
      dim = n_icu_stays
    ),
    beta_slope = array(
      data = rnorm(n = n_icu_stays * 2, mean = c(2500, 1500), sd = 500) %>% abs(),
      dim = c(n_icu_stays, 2)
    ),
    beta_zero = array(
      data = rnorm(n = n_icu_stays, mean = 7500, sd = 500),
      dim = n_icu_stays
    ),
    y_sigma = rnorm(n = 1, mean = 50, sd = 20) %>% abs()
  )
}

prefit <- stan_model('scripts/mimic-example/models/fluid-piecewise-linear.stan')
model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = 4,
  init = init_gen,
  control = list(
    max_treedepth = 12,
    adapt_delta = 0.99
  ),
  warmup = 1000,
  iter = 5e3
)

library(tidybayes)
library(bayesplot)

model_fit %>% 
  mcmc_trace(regex_pars = 'breakpoint_raw')

model_fit %>% 
  mcmc_combo(
    pars = c(
      'breakpoint_raw[3]',
      'breakpoint_raw[5]',
      'breakpoint_raw[9]',
      'breakpoint_raw[15]'
    ),
    combo = c('dens_overlay', 'trace')
  )

model_fit %>% 
  mcmc_trace(regex_pars = 'beta_zero')

model_fit %>% 
  mcmc_trace(regex_pars = 'beta_slope')

model_fit %>% 
  mcmc_trace(regex_pars = 'lp')

n_indiv <- stan_data$n_icu_stays

for (ii in 1 : n_indiv) {
  p1 <- model_fit %>% 
    mcmc_parcoord(
      pars = vars(
        param_range('breakpoint', ii),
        param_range('beta_zero', ii),
        param_glue('beta_slope[{indiv},{level}]', indiv = ii, level = 1 : 2),
        lp__
      ),
      transform = function(x) {(x - mean(x)) / sd(x)},
      np = nuts_params(model_fit),
      np_style = parcoord_style_np(div_size = 0.8, div_alpha = 0.8)
    ) +
    theme(
      plot.margin = margin(
        t = 0,
        r = 25,
        b = 10,
        l = 25
      )
    )
  
  plot_name <- sprintf('indiv-%02d.png', ii)
  
  ggsave(
    filename = file.path('plots', 'mimic-example', 'temp-parcoords', plot_name),
    plot = p1,
    width = 6,
    height = 8
  )
}

model_fit_long <- model_fit %>% 
  gather_draws(mu[n])

x_tbl <- tibble(
  n = 1 : stan_data$n_total_obs,
  x = res$x_obs,
  indiv_id = res$indiv_id
)

plot_tbl_by_chain <- model_fit_long %>% 
  group_by(n, .chain, .variable) %>% 
  point_interval(.point = mean, .exclude = c('.iteration', '.draw')) %>% 
  left_join(x_tbl, by = 'n') %>% 
  left_join(
    res %>% 
      select(indiv_id, plot_label),
    by = 'indiv_id'
  )

base_plot + 
  geom_ribbon(
    inherit.aes = FALSE,
    data = plot_tbl_by_chain %>% 
      mutate(.chain = as.factor(.chain)),
    mapping = aes(
      x = x,
      ymin = .lower,
      ymax = .upper,
      colour = .chain,
      group =
    ),
    alpha = 0.2
  )

