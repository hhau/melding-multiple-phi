library(tibble)
library(dplyr)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(data_seed)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

flog.info("surv-submodel-three: simulating data", name = base_filename)

simulated_data <- with(sim_settings, 
  bind_rows(lapply(1 : n_patients, function(patient_id) {
    beta_zero_true <- long_rand_ef[patient_id]

    obs_times <- sort(runif(
      n = n_obs_per_patient[patient_id]
    ))

    true_values <- beta_zero_true
    obs_values <- abs(
      true_values + 
      rnorm(n = n_obs_per_patient[patient_id], mean = 0, sd = 8 * sigma_noise)
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
