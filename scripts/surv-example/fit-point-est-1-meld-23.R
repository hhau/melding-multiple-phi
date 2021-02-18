library(magrittr)
library(rstan)
library(dplyr)
library(stringr)
library(parallel)
library(abind)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

flog.info("surv-fit-point-est-1-meld-23: reading data", name = base_filename)

# read in output from submodel 1 and 3, and submodel 2's data 
submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
submodel_three_output <- readRDS("rds/surv-example/submodel-three-output.rds")
submodel_two_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")
phi_12_post_median <- readRDS("rds/surv-example/stage-one-phi-12-posterior-median.rds")

submodel_one_samples <- submodel_one_output$samples
n_iter_submodel_one <- dim(submodel_one_samples)[1]
n_chain_submodel_one <- dim(submodel_one_samples)[2]

submodel_three_samples <- submodel_three_output$samples
n_iter_submodel_three <- dim(submodel_three_samples)[1]
n_chain_submodel_three <- dim(submodel_three_samples)[2]

stage_two_indices_names_psi_1 <- c(
  "submodel_one_iter_index", 
  "submodel_one_chain_index"
)

stage_two_indices_names_psi_3 <- c(
  "submodel_three_iter_index", 
  "submodel_three_chain_index"
)

# process into Stan data for all chains
phi_12_names <- c(
  sprintf("event_time[%d]", 1 : n_patients),
  sprintf("event_indicator[%d]", 1 : n_patients)
)

event_time_names <- grep("event_time", phi_12_names, value = TRUE)
event_indicator_names <- grep("event_indicator", phi_12_names, value = TRUE)

beta_gr <- expand.grid(1 : n_patients, 1 : n_long_beta)
phi_23_names <- sprintf("beta[%d,%d]", beta_gr[, 1], beta_gr[, 2])

n_phi_12 <- length(phi_12_names)
n_phi_23 <- length(phi_23_names)

stan_data <- list(
  n_patients = length(submodel_two_data$patient_id),
  n_long_beta = n_long_beta,
  baseline_measurement = submodel_two_data$baseline_val %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.numeric(),
  log_crude_event_rate = log(mean(submodel_one_samples[, , phi_12_names]))
)

flog.info("surv-fit-point-est-1-meld-23: compiling models", name = base_filename)

prefit_psi_two_step <- stan_model(
  "scripts/surv-example/models/submodel-two-psi-step.stan"
)

# NB: we can use the same stan file/fit for both phi_{1 \cap 2} and 
# phi_{2 \cap 3}, because we can control which parameters change / which
# are the same.
prefit_phi_step <- stan_model(
  "scripts/surv-example/models/submodel-two-phi-step.stan"
)

stanfit_phi_step <- sampling(
  prefit_phi_step,
  data = stan_data,
  chains = 1,
  cores = 1,
  iter = 1,
  refresh = 0
)

psi_two_names <- stanfit_phi_step@model_pars %>%
  grep(".*(theta|gamma|alpha).*", x = ., value = TRUE) %>% 
  magrittr::extract(!str_detect(., "long"))
  
n_psi_two_pars <- length(psi_two_names)

# general MCMC options
n_stage_two_chain <- 5
n_stage_two_iter <- 5e3 # 4e4 for min tail_ess in phi of ~500

list_res <- mclapply(1 : n_stage_two_chain, mc.cores = 5, function(chain_id) {
  # set up containers, remember we are abind'ing over the chains
  # phi - event times (no censoring yet)
  flog.info(
    sprintf("surv-fit-point-est-1-meld-23--chain-%d: initialising", chain_id), 
    name = base_filename
  )

  phi_23_samples <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, n_phi_23),
    dimnames = list(NULL, paste0("chain_", chain_id), phi_23_names)
  )
  
  psi_2_samples <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, n_psi_two_pars),
    dimnames = list(NULL, paste0("chain_", chain_id), psi_two_names)
  )

  psi_3_indices <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, length(stage_two_indices_names_psi_3)),
    dimnames = list(NULL, paste0("chain_", chain_id), stage_two_indices_names_psi_3)
  )
  
  # initialise containers
  psi_3_iter_index_init <- sample.int(n = n_iter_submodel_three, size = 1)
  psi_3_chain_index_init <- sample.int(n = n_chain_submodel_three, size = 1)
  psi_3_indices[1, 1, ] <- c(psi_3_iter_index_init, psi_3_chain_index_init)

  phi_23_samples[1, 1, ] <- submodel_three_samples[
    psi_3_iter_index_init,
    psi_3_chain_index_init, 
    phi_23_names
  ]
  
  psi_2_samples[1, 1, ] <- c(
    beta_zero = rnorm(n = 1, mean = 4, sd = 1),
    beta_one = rnorm(n = 1, mean = 0.5, sd = 0.2),
    hazard_gamma = abs(rnorm(n = 1, mean = 0.6, sd = 0.2)),
    alpha = rnorm(n = 1, mean = 0, sd = 0.1)
  )
  
  flog.info(
    sprintf("surv-fit-point-est-1-meld-23--chain-%d: initialised okay", chain_id), 
    name = base_filename
  )

  for (ii in 2 : n_stage_two_iter) {
    phi_12_list <- phi_12_post_median

    phi_23_curr_list <- list(
      long_beta = matrix(
        phi_23_samples[ii - 1, 1, ],
        ncol = n_long_beta,
        nrow = n_patients
      )
    )
    
    # psi step
    psi_step_data <- c(
      stan_data,
      phi_12_list,
      phi_23_curr_list
    )

    psi_step <- sampling(
      object = prefit_psi_two_step,
      data = psi_step_data,
      init = list(as.list(psi_2_samples[ii - 1, 1, ])),
      chains = 1,
      iter = 6,
      warmup = 5,
      refresh = 0
    )

    psi_2_samples[ii, 1, ] <- as.array(psi_step, pars = psi_two_names)
    psi_2_curr_list <- as.list(psi_2_samples[ii, 1, ])
    
    # phi_{2 \cap 3} step
    phi_23_iter_index_prop <- sample.int(n = n_iter_submodel_three, size = 1)
    phi_23_chain_index_prop <- sample.int(n = n_chain_submodel_three, size = 1)
    phi_23_prop_list <- list(
      long_beta = matrix(
        submodel_three_samples[
          phi_23_iter_index_prop, 
          phi_23_chain_index_prop, 
          phi_23_names
        ],
        ncol = n_long_beta,
        nrow = n_patients
      )
    )
    
    log_prob_phi_23_step_prop <- log_prob(
      object = stanfit_phi_step,
      upars = unconstrain_pars(
        stanfit_phi_step,
        pars = c(
          psi_2_curr_list,
          phi_12_list,
          phi_23_prop_list
        )
      )
    )

    log_prob_phi_23_step_curr <- log_prob(
      object = stanfit_phi_step,
      upars = unconstrain_pars(
        stanfit_phi_step,
        pars = c(
          psi_2_curr_list,
          phi_12_list,
          phi_23_curr_list
        )
      )
    )

    log_alpha_23_step <- log_prob_phi_23_step_prop - log_prob_phi_23_step_curr

    if (runif(1) < exp(log_alpha_23_step)) {
      psi_3_indices[ii, 1, ] <- c(phi_23_iter_index_prop, phi_23_chain_index_prop)
      phi_23_samples[ii, 1, ] <- phi_23_prop_list$long_beta
    } else {
      psi_3_indices[ii, 1, ] <- psi_3_indices[ii - 1, 1, ]
      phi_23_samples[ii, 1, ] <-  phi_23_samples[ii - 1, 1, ]
    }
    
    if (ii %% 500 == 0) {
      flog.info(
        sprintf("surv-fit-point-est-1-meld-23--chain-%d: iteration-%d", chain_id, ii), 
        name = base_filename
      )
    }
  }
  
  inner_res <- list(
    psi_2_samples = psi_2_samples,
    phi_23_samples = phi_23_samples,
    psi_3_indices = psi_3_indices
  )
  
  return(inner_res)
})

bind_named_sublists <- function(outer_list, name) {
  lapply(outer_list, function(sub_list) sub_list[[name]]) %>% 
    abind(along = 2)
}

psi_2_samples <- bind_named_sublists(list_res, "psi_2_samples")
phi_23_samples <- bind_named_sublists(list_res, "phi_23_samples")
psi_3_indices <- bind_named_sublists(list_res, "psi_3_indices")

flog.info("surv-fit-point-est-1-meld-23: writing to disk", name = base_filename)

saveRDS(
  object = phi_23_samples,
  file = "rds/surv-example/point-est-1-meld-23-phi-23-samples.rds"
)

saveRDS(
  object = psi_2_samples,
  file = "rds/surv-example/point-est-1-meld-23-psi-2-samples.rds"
)

saveRDS(
  object = psi_3_indices,
  file = "rds/surv-example/point-est-1-meld-23-psi-3-indices.rds"
)

flog.info("surv-fit-point-est-1-meld-23: finished writing", name = base_filename)
