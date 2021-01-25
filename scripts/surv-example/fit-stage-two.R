library(rstan)
library(dplyr)
library(parallel)
library(abind)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

flog.info("surv-fit-stage-two: reading data", name = base_filename)

# read in output from submodel 1 and submodel 2's data 
submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
submodel_two_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")

submodel_one_samples <- submodel_one_output$samples
n_iter_submodel_one <- dim(submodel_one_samples)[1]
n_chain_submodel_one <- dim(submodel_one_samples)[2]
stage_two_indices_names <- c(
  "submodel_one_iter_index", 
  "submodel_one_chain_index"
)

# process into Stan data for all chains
event_time_names <- sprintf("event_time[%d]", 1 : n_patients)
stan_data <- list(
  n_patients = length(submodel_two_data$patient_id),
  baseline_measurement = submodel_two_data$baseline_val %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.numeric(),
  log_crude_event_rate = log(mean(submodel_one_samples[, , event_time_names]))
)

flog.info("surv-fit-stage-two: compiling models", name = base_filename)

prefit_two_psi_step <- stan_model(
  "scripts/surv-example/models/submodel-two-psi-step.stan"
)

prefit_two_phi_step <- stan_model(
  "scripts/surv-example/models/submodel-two-phi-step.stan"
)

stanfit_two_phi_step <- sampling(
  prefit_two_phi_step,
  data = stan_data,
  chains = 1,
  cores = 1,
  iter = 1,
  refresh = 0
)

psi_two_names <- stanfit_two_phi_step@model_pars %>%
  grep(".*(beta|gamma).*", x = ., value = TRUE)

n_psi_two_pars <- length(psi_two_names)

# general MCMC options
n_stage_two_chain <- 5
n_stage_two_iter <- 4e3 # 4e4 for min tail_ess in phi of ~500

list_res <- mclapply(1 : n_stage_two_chain, mc.cores = 5, function(chain_id) {
  # set up containers, remember we are abind'ing over the chains
  # phi - event times (no censoring yet)
  flog.info(
    sprintf("surv-fit-stage-two--chain-%d: initialising", chain_id), 
    name = base_filename
  )

  phi_samples <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, n_patients),
    dimnames = list(NULL, paste0("chain_", chain_id), event_time_names)
  )
  
  psi_samples <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, n_psi_two_pars),
    dimnames = list(NULL, paste0("chain_", chain_id), psi_two_names)
  )
  
  indices <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, length(stage_two_indices_names)),
    dimnames = list(NULL, paste0("chain_", chain_id), stage_two_indices_names)
  )
  
  # initialise containers
  chain_index_init <- sample.int(n = n_chain_submodel_one, size = 1)
  iter_index_init <- sample.int(n = n_iter_submodel_one, size = 1)
  
  phi_samples[1, 1, ] <- submodel_one_samples[
    iter_index_init,
    chain_index_init, 
    event_time_names
  ]
  
  indices[1, 1, ] <- c(iter_index_init, chain_index_init)
  psi_samples[1, 1, ] <- c(
    beta_zero = rnorm(n = 1, mean = 4, sd = 1),
    beta_one = rnorm(n = 1, mean = 0.5, sd = 0.2),
    hazard_gamma = abs(rnorm(n = 1, mean = 0.6, sd = 0.2))
  )

  for (ii in 2 : n_stage_two_iter) {
    # psi step
    psi_step_data <- c(
      stan_data,
      event_time = list(as.numeric(phi_samples[ii - 1, 1, ]))
    )

    psi_step <- sampling(
      object = prefit_two_psi_step,
      data = psi_step_data,
      init = list(as.list(psi_samples[ii - 1, 1, ])),
      chains = 1,
      iter = 26,
      warmup = 25,
      refresh = 0
    )

    psi_samples[ii, 1, ] <- as.array(psi_step, pars = psi_two_names)
    psi_curr_list <- as.list(psi_samples[ii, 1, ])

    # phi step
    iter_index_prop <- sample.int(n = n_iter_submodel_one, size = 1)
    chain_index_prop <- sample.int(n = n_chain_submodel_one, size = 1)
    phi_prop_list <- list(event_time = as.numeric(
      submodel_one_samples[iter_index_prop, chain_index_prop, event_time_names]
    ))
    
    log_prob_phi_step_prop <- log_prob(
      object = stanfit_two_phi_step,
      upars = unconstrain_pars(
        stanfit_two_phi_step,
        pars = c(
          phi_prop_list,
          psi_curr_list
        )
      )
    )

    log_prob_phi_step_curr <- log_prob(
      object = stanfit_two_phi_step,
      upars = unconstrain_pars(
        stanfit_two_phi_step,
        pars = c(
          event_time = list(as.numeric(phi_samples[ii - 1, 1, ])),
          psi_curr_list
        )
      )
    )

    log_alpha <- log_prob_phi_step_prop - log_prob_phi_step_curr

    if (runif(1) < exp(log_alpha)) {
      indices[ii, 1, ] <- c(iter_index_prop, chain_index_prop)
      phi_samples[ii, 1, ] <- phi_prop_list$event_time
    } else {
      indices[ii, 1, ] <- indices[ii - 1, 1, ]
      phi_samples[ii, 1, ] <-  phi_samples[ii - 1, 1, ]
    }

    if (ii %% 500 == 0) {
      flog.info(
        sprintf("surv-fit-stage-two--chain-%d: iteration-%d", chain_id, ii), 
        name = base_filename
      )
    }
  }
  
  inner_res <- list(
    indices = indices,
    phi_samples = phi_samples,
    psi_samples = psi_samples
  )
})

bind_named_sublists <- function(outer_list, name) {
  lapply(outer_list, function(sub_list) sub_list[[name]]) %>% 
    abind(along = 2)
}

indices <- bind_named_sublists(list_res, "indices")
phi_samples <- bind_named_sublists(list_res, "phi_samples")
psi_samples <- bind_named_sublists(list_res, "psi_samples")

flog.info("surv-fit-stage-two: writing to disk", name = base_filename)

saveRDS(
  object = phi_samples,
  file = "rds/surv-example/stage-two-phi-samples.rds"
)

saveRDS(
  object = indices,
  file = "rds/surv-example/stage-two-indices.rds"
)

saveRDS(
  object = psi_samples,
  file = "rds/surv-example/stage-two-psi-samples.rds"
)

flog.info("surv-fit-stage-two: finished writing", name = base_filename)
