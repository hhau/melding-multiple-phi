source("scripts/surv-example/GLOBALS.R")

set.seed(data_seed)
# Submodel three settings.
n_obs_per_patient <- pmax(rpois(n_patients, 8), 3)
n_biomarkers <- 1

sigma_noise <- 0.1

# some plot settings
n_plot <- 100
x_plot <- seq(from = 0, to = 1, length.out = n_plot)

submodel_three_settings <- list(
  n_patients = n_patients,
  n_biomarkers = n_biomarkers, 
  n_obs_per_patient = n_obs_per_patient,
  sigma_noise = sigma_noise,
  n_plot = n_plot,
  x_plot = x_plot
)

saveRDS(
  object = submodel_three_settings,
  file = "rds/surv-example/submodel-three-simulation-settings.rds"
)
