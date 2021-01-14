library(tibble)
library(dplyr)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(data_seed)

# peek at submodel one's settings (non-generative process) to see which
# individuals have the event
submodel_one_settings <- readRDS(
  "rds/surv-example/submodel-one-simulation-settings.rds"
)

flog.info("surv-submodel-two: simulating data", name = base_filename)

simulated_data <- with(submodel_one_settings, {
  lapply(1 : n_patients, function(patient_id) {
    tibble(
      patient_id = patient_id,
      baseline_val = ifelse(
        event_indicator[patient_id],
        rnorm(n = 1, mean = 2, sd = 1),
        rnorm(n = 1, mean = 0, sd = 1)
      )
    )
  }) %>% bind_rows()
})

flog.info("surv-submodel-two: saving simulated data", name = base_filename)

saveRDS(
  file = "rds/surv-example/submodel-two-simulated-data.rds",
  object = simulated_data
)
