library(rstan)
library(bayesplot)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

# I should really copy mike and do this with argparse
simulated_data <- readRDS(
  file = "rds/surv-example/submodel-three-simulated-data.rds"
)

sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

flog.info("surv-fit-submodel-three: compiling model", name = base_filename)

model_prefit <- stan_model(
  file = "scripts/surv-example/models/submodel-three.stan"
)

stan_input_data <- with(simulated_data, 
  list(
    n_obs = length(patient_id),
    n_patients = length(unique(patient_id)),
    Y = measurement,
    obs_ids = patient_id,
    obs_times = time,
    n_plot = sim_settings$n_plot,
    x_plot = sim_settings$x_plot
  )
)

flog.info(
  "surv-fit-submodel-three: fitting model",
   name = base_filename
 )

model_fit <- sampling(
  model_prefit,
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

model_three_sampler_output <- list(
  samples = raw_samples,
  nuts_params = nuts_params_for_samples
)

flog.info(
  "surv-fit-submodel-three: saving output",
  name = base_filename
)

saveRDS(
  object = model_three_sampler_output,
  file = "rds/surv-example/submodel-three-output.rds"
)
