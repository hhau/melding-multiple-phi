library(tidybayes)
library(dplyr)

source("scripts/common/mcmc-util.R")

simulated_data <- readRDS(
  file = "rds/surv-example/submodel-one-simulated-data.rds"
)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

# read in model output
model_output <- readRDS(
  file = "rds/surv-example/submodel-one-output.rds"
)

submodel_one_point_tbl <- simulated_data %>% 
  mutate(
    model_name = "event_time",
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

submodel_one_fitted_tbl <- posterior_plot_data %>% 
  mutate(
    model_name = "event_time",
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

submodel_one_threshold_tbl <- submodel_one_fitted_tbl %>%
  filter(time %in% c(0, 1)) %>%
  mutate(
    y_mean = NA,
    y_lower = NA,
    y_upper = NA,
    threshold_value = sim_settings$y_threshold
  )

event_res <- model_output$samples %>%
  array_to_mcmc_list() %>%  
  spread_draws(event_time[patient_id])

event_data <- event_res %>% 
  point_interval(.width = 0.8)

trimmed_event_data <- event_data %>% 
  rowwise() %>% 
  mutate(
    across(.lower, function(x) ifelse(between(x, 0, 1), x, 0)),
    across(.upper, function(y) ifelse(between(y, 0, 1), y, 1))
  ) %>% 
  filter(
    !(.lower == 0 & .upper == 1)
  )

submodel_one_interval_tbl <- trimmed_event_data %>% 
  mutate(
    model_name = "event_time",
    time = NA,
    y_point = NA,
    y_mean = NA,
    y_lower = NA,
    y_upper = NA,
    interval_lower = .lower,
    interval_upper = .upper,
    threshold_value = NA
  ) %>%
  select(-c(event_time, .lower, .upper, .width, .point, .interval)) %>%
  left_join(
    submodel_one_point_tbl %>% select(patient_id, event_indicator),
    by = c("patient_id")
  )

submodel_one_tbl <- bind_rows(
  submodel_one_point_tbl,
  submodel_one_fitted_tbl,
  submodel_one_threshold_tbl,
  submodel_one_interval_tbl
) %>%
  mutate(event_indicator = as.factor(event_indicator))

saveRDS(
  object = submodel_one_tbl,
  file = "rds/surv-example/plot-submodel-1-data.rds"
)
