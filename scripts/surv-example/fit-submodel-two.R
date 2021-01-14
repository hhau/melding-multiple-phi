# technically this isn't fitting submodel two as a standalone submodel.
# I could do that, but that would be a very standard survival model, with a
# fixed event time.
# this is more like what I want to do

# libs and logs
library(parallel)
library(rstan)
library(bayesplot)
library(rstanarm)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

# read in output from submodel 1 and submodel 2's data 
submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
submodel_two_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")

# compile the necessary models and generate stanfit object
psi_step_stan_model <- stan_model("scripts/surv-example/models/submodel-two-psi-step.stan")

# general MCMC options

# second stage MCMC, note that because we aren't using, from the first stage, 
# something that would typically be called a parameter, we can actually just
# let Stan do all of the fitting. this also means we could do
## no this is all wrong
## The above would be some kind of cut-like posterior,
## but that's not what we're doing. 
## the baseline covariates / parameter will inform which event times
## were plausible. Whether or not this makes any sense (physically) is 
## debatable, but it is what _should_ happen here given the melding process.


## some temporary code to see if this has a hope of working
event_time_names <- sprintf("event_time[%d]", 1 : n_patients)

stan_data <- list(
  n_patients = length(submodel_two_data$patient_id),
  baseline_measurement = submodel_two_data$baseline_val,
  event_time = submodel_one_output$samples[
    sample(1 : dim(submodel_one_output$samples)[1], 1), 
    sample(1 : dim(submodel_one_output$samples)[2], 1), 
    event_time_names
  ]
)

surv_reg_df <- data.frame(
  event_time = stan_data$event_time,
  event_indicator = rep(1, stan_data$n_patients),
  baseline = stan_data$baseline_measurement
)

test_fit <- sampling(
  psi_step_stan_model,
  data = stan_data,
  cores = 4
)

res_surv_stan <- stan_surv(
  Surv(event_time, event_indicator) ~ baseline, 
  data = surv_reg_df, 
  basehaz = "weibull",
  cores = 4,
  refresh = 0
)

res_surv_stan

test_samples <- rstan::extract(test_fit, permuted = FALSE, inc_warmup = TRUE)
plot_pars <- c("lp__", "beta_zero", "beta_one", "hazard_lambda")
mcmc_trace(test_samples[, , plot_pars])


sink(file = "rstanarm_surv_stan_code.stan")
res$stanfit@stanmodel
sink()
