library(tibble)
library(dplyr)

source("scripts/surv-example/GLOBALS.R")

set.seed(data_seed)

# Submodel one settings.
event_prop <- 0.5
n_biomarkers <- 1

n_obs_per_patient <- sample(2 : 8, size = n_patients, replace = TRUE)
event_indicator <- rbinom(n_patients, 1, event_prop)

sigma_latent <- 1
sigma_noise <- 0.1
y_threshold <- 0.2

# some plot settings
n_plot <- 100
x_plot <- seq(from = 0, to = 1, length.out = n_plot)

n_latent_params <- 2 + n_long_beta

no_event_cov_mat <- diag(n_latent_params)
event_cov_mat <- matrix(
  data = rep(0.5, times = n_latent_params^2), 
  ncol = n_latent_params,
  nrow = n_latent_params
) + 
  diag(sigma_latent - 0.5, n_latent_params)

event_cov_chol <- chol(event_cov_mat)

global_data <- lapply(1 : n_patients, function(patient_id) {
  noise_obs <- matrix(
    rnorm(n = n_latent_params), 
    ncol = n_latent_params, 
    nrow = 1
  )

  if (event_indicator[patient_id] == 1) {
    vals <- c(3, 5, 3, 1.5) + noise_obs %*% event_cov_chol
  } else {
    vals <- noise_obs
  }

  res <- tibble(
    event_model_slope = -abs(vals[1, 1]) / 2,
    baseline_cov = vals[1, 2],
    long_rand_ef = vals[1, 3],
    long_rand_slope = vals[1, 4]
  )
}) %>%
  bind_rows()


settings_and_joint_data <- list(
  n_patients = n_patients, 
  event_prop = event_prop,
  n_obs_per_patient = n_obs_per_patient,
  event_model_slope = global_data$event_model_slope,
  baseline_cov = global_data$baseline_cov,
  long_rand_ef = global_data$long_rand_ef,
  long_rand_slope = global_data$long_rand_slope,
  event_indicator = event_indicator,
  sigma_noise = sigma_noise,
  y_threshold = y_threshold,
  n_plot = n_plot,
  x_plot = x_plot
)

saveRDS(
  object = settings_and_joint_data,
  file = "rds/surv-example/simulation-settings-and-joint-data.rds"
)
