library(rstan)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")
set.seed(sim_seed)

flog.info("Surv ex - point est approx: loading data", name = base_filename)

submodel_one_samples <- readRDS("rds/surv-example/stage-one-phi-12-samples.rds")
submodel_two_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")
phi_12_post_median <- readRDS("rds/surv-example/stage-one-phi-12-posterior-median.rds")
phi_23_post_median <- readRDS("rds/surv-example/stage-one-phi-23-posterior-median.rds")

flog.info("Surv ex - point est approx: compiling", name = base_filename)

prefit <- stan_model(
  "scripts/surv-example/models/submodel-two-psi-step.stan"
)

stan_data <- list(
  n_patients = length(submodel_two_data$patient_id),
  n_long_beta = length(phi_23_post_median), # there should be a list element per beta
  baseline_measurement = submodel_two_data$baseline_val %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.numeric(),
  log_crude_event_rate = log(mean(submodel_one_samples)), # annoying but we should be consistent
  event_indicator = phi_12_post_median$event_indicator,
  event_time = phi_12_post_median$event_time,
  long_beta = matrix(
    phi_23_post_median$long_beta,
    ncol = length(phi_23_post_median)
  )
)

flog.info("Surv ex - point est approx: sampling", name = base_filename)

model_fit <- sampling(
  prefit,
  data = stan_data,
  include = FALSE,
  pars = "z_common",
  # cores = 5,
  chains = 5,
  iter = 6e3,
  warmup = 1e3,
  control = list(adapt_delta = 0.95),
  refresh = 0
)

samples <- as.array(model_fit)

flog.info("Surv ex - point est approx: saving samples", name = base_filename)

saveRDS(
  file = "rds/surv-example/point-est-psi-2-samples.rds",
  object = samples
)
