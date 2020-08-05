library(tibble)
library(dplyr)

source("scripts/common/logger-setup.R")

sim_settings <- readRDS("rds/surv-example/submodel-one-simulation-settings.rds")

# TODO: censoring time / censoring in general?
# will this lead to truncated distributions?
# Possibly too hard

flog.info("surv: simulating data")

simulated_data <- with(sim_settings, 
  bind_rows(lapply(1 : n_patients, function(patient_id) {
    alpha_zero_true <- rnorm(n = 1, mean = 1, sd = 0.1)
    alpha_one_true <- rnorm(n = 1, mean = 0.5, sd = 0.1)
    beta_one_true <- rnorm(
      n = 1,
      mean = ifelse(event_indicator[patient_id], -1, 0),
      sd = ifelse(event_indicator[patient_id], 0.2, 0.1)
    )
    obs_times <- sort(runif(n = n_obs_per_patient[patient_id]))
    covariate_values <- rnorm(n = n_obs_per_patient[patient_id], sd = 0.25)

    true_values <- alpha_zero_true +
     alpha_one_true * covariate_values + 
     beta_one_true * obs_times

    regression_errors <- rnorm(
      n = n_obs_per_patient[patient_id],
      mean = 0,
      sd = sigma_noise_x
    )

    trajectory_errors <- rnorm(
      n = n_obs_per_patient[patient_id],
      mean = 0,
      sd = sigma_noise_t
    )

    obs_values <- true_values + regression_errors + trajectory_errors

    res <- tibble(
      patient_id = patient_id,
      time = obs_times,
      measurement = obs_values,
      covariate_values = covariate_values,
      event_indicator = event_indicator[patient_id]
    )
  }))
)

flog.info("surv: saving simulated data to disk")

saveRDS(
  object = simulated_data,
  file = "rds/surv-example/submodel-one-simulated-data.rds"
)
