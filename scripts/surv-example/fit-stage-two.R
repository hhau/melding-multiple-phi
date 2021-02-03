library(magrittr)
library(rstan)
library(dplyr)
library(stringr)
library(parallel)
library(abind)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

flog.info("surv-fit-stage-two: reading data", name = base_filename)

# read in output from submodel 1 and 3, and submodel 2's data 
submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
submodel_three_output <- readRDS("rds/surv-example/submodel-three-output.rds")
submodel_two_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")

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
phi_12_names <- sprintf("event_time[%d]", 1 : n_patients)
phi_23_names <- c(
  sprintf("beta_zero[%d]", 1 : n_patients)#,
  # sprintf("beta_one[%d]", 1 : n_patients)
)

n_phi_12 <- length(phi_12_names)
n_phi_23 <- length(phi_23_names)
n_long_beta <- 1 # 1 = intercept only model, 2 = intercept + slope

stan_data <- list(
  n_patients = length(submodel_two_data$patient_id),
  n_long_beta = n_long_beta,
  baseline_measurement = submodel_two_data$baseline_val %>%
    scale(center = TRUE, scale = FALSE) %>%
    as.numeric(),
  log_crude_event_rate = log(mean(submodel_one_samples[, , phi_12_names]))
)

flog.info("surv-fit-stage-two: compiling models", name = base_filename)

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
  grep(".*(beta|gamma|alpha).*", x = ., value = TRUE) %>% 
  magrittr::extract(!str_detect(., "long"))
  
n_psi_two_pars <- length(psi_two_names)

# general MCMC options
n_stage_two_chain <- 5
n_stage_two_iter <- 5e3 # 4e4 for min tail_ess in phi of ~500

list_res <- mclapply(1 : n_stage_two_chain, mc.cores = 5, function(chain_id) {
  # set up containers, remember we are abind'ing over the chains
  # phi - event times (no censoring yet)
  flog.info(
    sprintf("surv-fit-stage-two--chain-%d: initialising", chain_id), 
    name = base_filename
  )

  phi_12_samples <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, n_phi_12),
    dimnames = list(NULL, paste0("chain_", chain_id), phi_12_names)
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
  
  psi_1_indices <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, length(stage_two_indices_names_psi_1)),
    dimnames = list(NULL, paste0("chain_", chain_id), stage_two_indices_names_psi_1)
  )

  psi_3_indices <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, length(stage_two_indices_names_psi_3)),
    dimnames = list(NULL, paste0("chain_", chain_id), stage_two_indices_names_psi_3)
  )
  
  # initialise containers
  psi_1_iter_index_init <- sample.int(n = n_iter_submodel_one, size = 1)
  psi_1_chain_index_init <- sample.int(n = n_chain_submodel_one, size = 1)
  psi_1_indices[1, 1, ] <- c(psi_1_iter_index_init, psi_1_chain_index_init)

  psi_3_iter_index_init <- sample.int(n = n_iter_submodel_three, size = 1)
  psi_3_chain_index_init <- sample.int(n = n_chain_submodel_three, size = 1)
  psi_3_indices[1, 1, ] <- c(psi_3_iter_index_init, psi_3_chain_index_init)
  
  phi_12_samples[1, 1, ] <- submodel_one_samples[
    psi_1_iter_index_init,
    psi_1_chain_index_init, 
    phi_12_names
  ]

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
    sprintf("surv-fit-stage-two--chain-%d: initialised okay", chain_id), 
    name = base_filename
  )

  for (ii in 2 : n_stage_two_iter) {
    phi_12_curr_list <- list(as.numeric(phi_12_samples[ii - 1, 1, ]))
    phi_23_curr_list <- list(as.numeric(phi_23_samples[ii - 1, 1, ]))
    
    # psi step
    psi_step_data <- c(
      stan_data,
      event_time = list(as.numeric(phi_12_samples[ii - 1, 1, ])),
      long_beta = list(
        matrix(
          data = phi_23_samples[ii - 1, 1, ],
          nrow = stan_data$n_patients,
          ncol = n_long_beta
        )
      )
    )

    psi_step <- sampling(
      object = prefit_psi_two_step,
      data = psi_step_data,
      # init = list(as.list(psi_2_samples[ii - 1, 1, ])),
      chains = 1,
      iter = 21,
      warmup = 20,
      refresh = 0
    )

    psi_2_samples[ii, 1, ] <- as.array(psi_step, pars = psi_two_names)
    psi_2_curr_list <- as.list(psi_2_samples[ii, 1, ])

    # phi_{1 \cap 2} step
    phi_12_iter_index_prop <- sample.int(n = n_iter_submodel_one, size = 1)
    phi_12_chain_index_prop <- sample.int(n = n_chain_submodel_one, size = 1)
    phi_12_prop_list <- list(event_time = as.numeric(
      submodel_one_samples[
        phi_12_iter_index_prop, 
        phi_12_chain_index_prop, 
        phi_12_names
      ]
    ))
    
    log_prob_phi_12_step_prop <- log_prob(
      object = stanfit_phi_step,
      upars = unconstrain_pars(
        stanfit_phi_step,
        pars = c(
          psi_2_curr_list,
          phi_12_prop_list,
          long_beta = phi_23_curr_list
        )
      )
    )

    log_prob_phi_12_step_curr <- log_prob(
      object = stanfit_phi_step,
      upars = unconstrain_pars(
        stanfit_phi_step,
        pars = c(
          psi_2_curr_list,
          event_time = phi_12_curr_list,
          long_beta = phi_23_curr_list
        )
      )
    )

    log_alpha_12_step <- log_prob_phi_12_step_prop - log_prob_phi_12_step_curr

    if (runif(1) < exp(log_alpha_12_step)) {
      psi_1_indices[ii, 1, ] <- c(phi_12_iter_index_prop, phi_12_chain_index_prop)
      phi_12_samples[ii, 1, ] <- phi_12_prop_list$event_time
      phi_12_curr_list <- list(phi_12_prop_list$event_time) # for phi_{2 \cap 3} step
    } else {
      psi_1_indices[ii, 1, ] <- psi_1_indices[ii - 1, 1, ]
      phi_12_samples[ii, 1, ] <-  phi_12_samples[ii - 1, 1, ]
    }
    
    # phi_{2 \cap 3} step
    phi_23_iter_index_prop <- sample.int(n = n_iter_submodel_three, size = 1)
    phi_23_chain_index_prop <- sample.int(n = n_chain_submodel_three, size = 1)
    phi_23_prop_list <- list(long_beta = as.numeric(
      submodel_three_samples[
        phi_23_iter_index_prop, 
        phi_23_chain_index_prop, 
        phi_23_names
      ]
    ))
    
    log_prob_phi_23_step_prop <- log_prob(
      object = stanfit_phi_step,
      upars = unconstrain_pars(
        stanfit_phi_step,
        pars = c(
          psi_2_curr_list,
          event_time = phi_12_curr_list,
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
          event_time = phi_12_curr_list,
          long_beta = phi_23_curr_list
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
        sprintf("surv-fit-stage-two--chain-%d: iteration-%d", chain_id, ii), 
        name = base_filename
      )
    }
  }
  
  inner_res <- list(
    psi_1_indices = psi_1_indices,
    phi_12_samples = phi_12_samples,
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

psi_1_indices <- bind_named_sublists(list_res, "psi_1_indices")
phi_12_samples <- bind_named_sublists(list_res, "phi_12_samples")
psi_2_samples <- bind_named_sublists(list_res, "psi_2_samples")
phi_23_samples <- bind_named_sublists(list_res, "phi_23_samples")
psi_3_indices <- bind_named_sublists(list_res, "psi_3_indices")

flog.info("surv-fit-stage-two: writing to disk", name = base_filename)

saveRDS(
  object = phi_12_samples,
  file = "rds/surv-example/stage-two-phi-12-samples.rds"
)

saveRDS(
  object = phi_23_samples,
  file = "rds/surv-example/stage-two-phi-23-samples.rds"
)

saveRDS(
  object = psi_2_samples,
  file = "rds/surv-example/stage-two-psi-2-samples.rds"
)

saveRDS(
  object = psi_1_indices,
  file = "rds/surv-example/stage-two-psi-1-indices.rds"
)

saveRDS(
  object = psi_3_indices,
  file = "rds/surv-example/stage-two-psi-3-indices.rds"
)


flog.info("surv-fit-stage-two: finished writing", name = base_filename)
