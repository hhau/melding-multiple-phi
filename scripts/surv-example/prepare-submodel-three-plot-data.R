library(tidybayes)
library(dplyr)

source("scripts/common/mcmc-util.R")

# read in data
simulated_data <- readRDS(
  file = "rds/surv-example/submodel-three-simulated-data.rds"
)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

# read in model output
model_output <- readRDS(
  file = "rds/surv-example/submodel-three-output.rds"
)

submodel_three_point_tbl <- simulated_data %>% 
  mutate(
    model_name = "longitudinal",
    y_point = measurement,
    y_mean = NA,
    y_lower = NULL,
    y_upper = NA,
    interval_lower = NA,
    interval_upper = NA,
    threshold_value = NA
  ) %>% 
  select(-c(measurement))

res <- model_output$samples %>%
  array_to_mcmc_list() %>%  
  spread_draws(plot_mu[patient_id, plot_x])

res$plot_x <- sim_settings$x_plot[res$plot_x]

posterior_plot_data <- res %>% 
  point_interval(
    .width = 0.8,
    .point = mean,
    .interval = qi
  )

event_df <- tibble(
  patient_id = 1 : sim_settings$n_patients,
  event_indicator = sim_settings$event_indicator
)

posterior_plot_data <- posterior_plot_data %>%
  left_join(event_df, by = c("patient_id" = "patient_id"))

submodel_three_fitted_tbl <- posterior_plot_data %>% 
  mutate(
    model_name = "longitudinal",
    time = plot_x,
    y_point = NA,
    y_mean = plot_mu,
    y_lower = .lower,
    y_upper = .upper,
    interval_lower = NA,
    interval_upper = NA,
    threshold_value = NA
  ) %>% 
  select(-c(.width, .interval, plot_x, plot_mu, .lower, .upper, .point))

submodel_three_tbl <- bind_rows(
  submodel_three_point_tbl,
  submodel_three_fitted_tbl
) %>%
  mutate(event_indicator = as.factor(event_indicator))

saveRDS(
  object = submodel_three_tbl,
  file = "rds/surv-example/plot-submodel-3-data.rds"
)