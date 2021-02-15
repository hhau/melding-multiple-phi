library(rstan)
library(bayesplot)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

# I should really copy mike and do this with argparse
simulated_data <- readRDS(
  file = "rds/surv-example/submodel-one-simulated-data.rds"
)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

flog.info(
  "surv-fit-submodel-one: compiling model", 
  name = base_filename
)

model_prefit <- stan_model(
  file = "scripts/surv-example/models/submodel-one.stan"
)

stan_input_data <- with(simulated_data, 
  list(
    n_obs = length(patient_id),
    n_patients = length(unique(patient_id)),
    Y = measurement,
    obs_ids = patient_id,
    obs_times = time,
    y_threshold = sim_settings$y_threshold,
    n_plot = sim_settings$n_plot,
    x_plot = sim_settings$x_plot
  )
)

flog.info(
  "surv-fit-submodel-one: fitting model",
   name = base_filename
 )

model_fit <- sampling(
  model_prefit,
  seed = data_seed,
  data = stan_input_data,
  chains = 5,
  iter = 6e3,
  warmup = 1e3,
  cores = 5,
  control = list(
    adapt_delta = 0.95
  )
)

raw_samples <- as.array(model_fit)
nuts_params_for_samples <- nuts_params(model_fit)

model_one_sampler_output <- list(
  samples = raw_samples,
  nuts_params = nuts_params_for_samples
)

flog.info(
  "surv-fit-submodel-one: saving output",
  name = base_filename
)

saveRDS(
  object = model_one_sampler_output,
  file = "rds/surv-example/submodel-one-output.rds"
)
