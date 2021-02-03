library(tibble)
library(dplyr)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(data_seed)

sim_settings <- readRDS("rds/surv-example/submodel-three-simulation-settings.rds")
submodel_one_settings <- readRDS(
  "rds/surv-example/submodel-one-simulation-settings.rds"
)

flog.info("surv-submodel-three: simulating data", name = base_filename)

simulated_data <- with(sim_settings, 
  bind_rows(lapply(1 : n_patients, function(patient_id) {
    beta_zero_true <- ifelse(
      submodel_one_settings$event_indicator[patient_id],
      rnorm(n = 1, mean = 1.5, sd = 1),
      rnorm(n = 1, mean = 0, sd = 0.5)
    )

    obs_times <- sort(runif(
      n = n_obs_per_patient[patient_id]
    ))

    true_values <- beta_zero_true
    obs_values <- abs(
      true_values + 
      rnorm(n = n_obs_per_patient[patient_id], mean = 0, sd = sigma_noise)
    )
    res <- tibble(
      patient_id = patient_id,
      time = obs_times,
      measurement = obs_values
    )
  }))
)

flog.info(
  "surv-submodel-three: saving simulated data to disk",
  name = base_filename
)

saveRDS(
  object = simulated_data,
  file = "rds/surv-example/submodel-three-simulated-data.rds"
)
