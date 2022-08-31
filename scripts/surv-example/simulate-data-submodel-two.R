library(tibble)
library(dplyr)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(data_seed)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

flog.info("surv-submodel-two: simulating data", name = base_filename)

simulated_data <- with(sim_settings, {
  lapply(1 : n_patients, function(patient_id) {
    tibble(
      patient_id = patient_id,
      baseline_val = baseline_cov[patient_id]
    )
  }) %>% bind_rows()
})

flog.info("surv-submodel-two: saving simulated data", name = base_filename)

saveRDS(
  file = "rds/surv-example/submodel-two-simulated-data.rds",
  object = simulated_data
)
