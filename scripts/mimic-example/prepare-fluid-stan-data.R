library(dplyr)
library(purrr)
library(magrittr)

summarised_data <- readRDS("rds/mimic-example/combined-pf-and-summarised-fluids.rds")

cumulative_fluid_data <- summarised_data %>%
  filter(value_type == 'fluids') %>%
  group_by(icustay_id) %>%
  arrange(icustay_id, time_since_icu_adm) %>%
  mutate(
    cumulative_value = cumsum(value)
  )

n_icu_stays <- cumulative_fluid_data %>% 
  pull(icustay_id) %>% 
  n_distinct()

n_total_obs <- cumulative_fluid_data %>% 
  nrow()

subset_vector <- cumulative_fluid_data %>% 
  ungroup() %>% 
  group_split(icustay_id) %>% 
  map_int(~ nrow(.)) %>% 
  cumsum() %>% 
  add(1) %>% 
  c(1, .)

breakpoint_eps <- 0.01
breakpoint_lower_limits <- cumulative_fluid_data %>% 
  group_by(icustay_id) %>%
  summarise(lower = min(time_since_icu_adm)) %>% 
  pull(lower)

breakpoint_upper_limits <- cumulative_fluid_data %>% 
  group_by(icustay_id) %>% 
  summarise(upper = max(time_since_icu_adm)) %>% 
  pull(upper)

# breakpoint_ranges <- breakpoint_upper_limits - breakpoint_lower_limits
# breakpoint_lower_limits <- breakpoint_lower_limits + breakpoint_ranges * breakpoint_eps
# breakpoint_upper_limits <- breakpoint_upper_limits - breakpoint_ranges * breakpoint_eps

emp_intercept_stats <- cumulative_fluid_data %>% 
  group_by(icustay_id) %>% 
  summarise(
    mean = mean(cumulative_value),
    sd = sd(cumulative_value)
  )

stan_data <- list(
  n_icu_stays = n_icu_stays,
  n_total_obs = n_total_obs,
  subset_vector = subset_vector,
  breakpoint_lower = breakpoint_lower_limits,
  breakpoint_upper = breakpoint_upper_limits,
  y_vec = cumulative_fluid_data$cumulative_value,
  x_vec = cumulative_fluid_data$time_since_icu_adm,
  intercept_means = emp_intercept_stats$mean,
  intercept_sds = emp_intercept_stats$sd
)

saveRDS(
  file = 'rds/mimic-example/submodel-3-fluid-data-stan.rds',
  object = stan_data
)

library(rstan)

prefit <- stan_model('scripts/mimic-example/models/fluid-piecewise-linear.stan')
model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = 4,
  control = list(
    adapt_delta = 0.8,
    max_treedepth = 8
  )
)  
  
library(bayesplot)
library(tidybayes)

model_fit %>% 
  mcmc_trace(regex_pars = 'breakpoint')

model_fit %>% 
  mcmc_trace(regex_pars = 'beta_zero')

model_fit %>% 
  mcmc_trace(regex_pars = 'beta_slope')

model_fit %>% 
  mcmc_trace(regex_pars = 'lp')


model_fit_long <- model_fit %>% 
  gather_draws(mu[n])

x_tbl <- tibble(
  n = 1 : n_total_obs,
  x = cumulative_fluid_data$time_since_icu_adm,
  icustay_id = cumulative_fluid_data$icustay_id
)

plot_tbl_by_chain <- model_fit_long %>% 
  group_by(n, .chain, .variable) %>% 
  point_interval(.point = mean, .exclude = c('.iteration', '.draw')) %>% 
  left_join(x_tbl, by = 'n')

library(ggplot2)

base_plot <- ggplot(cumulative_fluid_data, aes(x = time_since_icu_adm, y = cumulative_value)) +
  geom_point() +
  facet_wrap(vars(icustay_id))

base_plot +
  geom_ribbon(
    inherit.aes = FALSE,
    data = plot_tbl_by_chain %>% 
      mutate(.chain = as.factor(.chain)),
    mapping = aes(
      x = x,
      ymin = .lower,
      ymax = .upper,
      colour = .chain
    ),
    alpha = 0.2
  )
