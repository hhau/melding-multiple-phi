source("scripts/surv-example/GLOBALS.R")

# Submodel one settings.
event_prop <- 0.5

n_obs_per_patient <- pmax(rpois(n_patients, 8), 3)
event_indicator <- rbinom(n_patients, 1, event_prop)

sigma_noise <- 0.1
y_threshold <- 0.2

# some plot settings
n_plot <- 100
x_plot <- seq(from = 0, to = 1, length.out = n_plot)

submodel_one_settings <- list(
  n_patients = n_patients, 
  event_prop = event_prop,
  n_obs_per_patient = n_obs_per_patient,
  event_indicator = event_indicator,
  sigma_noise = sigma_noise,
  y_threshold = y_threshold,
  n_plot = n_plot,
  x_plot = x_plot
)

saveRDS(
  object = submodel_one_settings,
  file = "rds/surv-example/submodel-one-simulation-settings.rds"
)
