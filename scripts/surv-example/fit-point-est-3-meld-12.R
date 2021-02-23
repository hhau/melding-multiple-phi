library(magrittr)
library(rstan)
library(dplyr)
library(stringr)
library(parallel)
library(abind)

source("scripts/common/logger-setup.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

flog.info("surv-fit-point-est-3-meld-12: reading data", name = base_filename)

# read in output from submodel 1 and 3, and submodel 2's data 
submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
submodel_two_data <- readRDS("rds/surv-example/submodel-two-simulated-data.rds")
phi_23_post_median <- readRDS("rds/surv-example/stage-one-phi-23-posterior-median.rds")
phi_23_post_median_mat <- matrix(
  unlist(phi_23_post_median),
  nrow = n_patients,
  ncol = n_long_beta,
  dimnames = list(
    NULL,
    names(phi_23_post_median)
  )
)

submodel_one_samples <- submodel_one_output$samples
n_iter_submodel_one <- dim(submodel_one_samples)[1]
n_chain_submodel_one <- dim(submodel_one_samples)[2]

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

flog.info("surv-fit-point-est-3-meld-12: compiling models", name = base_filename)

prefit_psi_two_step <- stan_model(
  "scripts/surv-example/models/submodel-two-psi-step.stan"
)

prefit_phi_indiv_step <- stan_model(
  "scripts/surv-example/models/submodel-two-phi-step-indiv.stan"
)

stanfit_phi_indiv_step <- sampling(
  prefit_phi_indiv_step,
  chains = 1,
  cores = 1,
  iter = 1,
  refresh = 0
)

psi_two_names <- c("theta_zero", "theta_one", "hazard_gamma", "alpha")
n_psi_two_pars <- length(psi_two_names)

# general MCMC options
n_stage_two_chain <- 5
n_stage_two_iter <- 5e3 # 4e4 for min tail_ess in phi of ~500

list_res <- mclapply(1 : n_stage_two_chain, mc.cores = 5, function(chain_id) {
  # set up containers, remember we are abind'ing over the chains
  # phi - event times (no censoring yet)
  flog.info(
    sprintf("surv-fit-point-est-3-meld-12--chain-%d: initialising", chain_id), 
    name = base_filename
  )

  phi_12_samples <- array(
    data = NA,
    dim = c(n_stage_two_iter, 1, n_phi_12),
    dimnames = list(NULL, paste0("chain_", chain_id), phi_12_names)
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
  
  # initialise containers
  psi_1_iter_index_init <- sample.int(n = n_iter_submodel_one, size = 1)
  psi_1_chain_index_init <- sample.int(n = n_chain_submodel_one, size = 1)
  psi_1_indices[1, 1, ] <- c(psi_1_iter_index_init, psi_1_chain_index_init)
  
  phi_12_samples[1, 1, ] <- submodel_one_samples[
    psi_1_iter_index_init,
    psi_1_chain_index_init, 
    phi_12_names
  ]
  
  psi_2_samples[1, 1, ] <- c(
    beta_zero = rnorm(n = 1, mean = 4, sd = 1),
    beta_one = rnorm(n = 1, mean = 0.5, sd = 0.2),
    hazard_gamma = abs(rnorm(n = 1, mean = 0.6, sd = 0.2)),
    alpha = rnorm(n = 1, mean = 0, sd = 0.1)
  )
  
  flog.info(
    sprintf("surv-fit-point-est-3-meld-12--chain-%d: initialised okay", chain_id), 
    name = base_filename
  )
  
  phi_23_list <- list(
    long_beta = matrix(
      data = c(
        phi_23_post_median$long_beta_zero, 
        phi_23_post_median$long_beta_one
      ),
      ncol = n_long_beta,
      nrow = n_patients
    )
  )

  for (ii in 2 : n_stage_two_iter) {
    phi_12_curr_list <- list(
      event_time = as.numeric(
        phi_12_samples[ii - 1, 1, event_time_names]
      ),
      event_indicator = as.integer(
        phi_12_samples[ii - 1, 1, event_indicator_names]
      )
    )
    
    # psi step
    psi_step_data <- c(
      stan_data,
      phi_12_curr_list,
      phi_23_list 
    )

    psi_step <- sampling(
      object = prefit_psi_two_step,
      data = psi_step_data,
      init = list(as.list(psi_2_samples[ii - 1, 1, ])),
      chains = 1,
      include = TRUE,
      pars = psi_two_names,
      iter = 6,
      warmup = 5,
      refresh = 0
    )

    psi_2_samples[ii, 1, ] <- as.array(psi_step, pars = psi_two_names)
    psi_2_curr_list <- as.list(psi_2_samples[ii, 1, ])

    # phi_{1 \cap 2} step
    # outside, randomly sample a scan order 
    # function / abstract out the updating? even more assigns / accesses
    scan_order <- sample.int(n_patients, n_patients, replace = FALSE)
    for (jj in 1 : n_patients) {
      indiv_index <- scan_order[jj]
      phi_12_iter_index_prop <- sample.int(n = n_iter_submodel_one, size = 1)
      phi_12_chain_index_prop <- sample.int(n = n_chain_submodel_one, size = 1)
   
      # define the names
      phi_12_indiv_names <- phi_12_names[c(indiv_index, indiv_index + n_patients)]
      
      # pack/denote the parameters
      phi_12_indiv_prop_list <- list(
        event_time = as.numeric(
          submodel_one_samples[
            phi_12_iter_index_prop,
            phi_12_chain_index_prop,
            phi_12_indiv_names[1]
          ]
        ),
        event_indicator = as.integer(
          submodel_one_samples[
            phi_12_iter_index_prop,
            phi_12_chain_index_prop,
            phi_12_indiv_names[2]
          ]
        ),
        eta_one = phi_23_post_median_mat[indiv_index, 2]
      )
      
      phi_12_indiv_curr_list <- list(
        event_time = phi_12_samples[ii - 1, 1, phi_12_indiv_names[1]],
        event_indicator = phi_12_samples[ii - 1, 1, phi_12_indiv_names[2]],
        eta_one = phi_23_post_median_mat[indiv_index, 2]
      )
      # eval the log probs
      log_prob_phi_12_indiv_prop <- log_prob(
        stanfit_phi_indiv_step,
        upars = unconstrain_pars(
          stanfit_phi_indiv_step,
          pars = c(psi_2_curr_list, phi_12_indiv_prop_list)
        )
      )
      
      log_prob_phi_12_indiv_curr <- log_prob(
        stanfit_phi_indiv_step,
        upars = unconstrain_pars(
          stanfit_phi_indiv_step,
          pars = c(psi_2_curr_list, phi_12_indiv_curr_list)
        )
      )
      
      log_alpha_phi_12_indiv <- log_prob_phi_12_indiv_prop - log_prob_phi_12_indiv_curr
      
      if (runif(1) < exp(log_alpha_phi_12_indiv)) {
        if (jj == 1) {
          psi_1_indices[ii, 1, ] <- c(phi_12_iter_index_prop, phi_12_chain_index_prop)
        }
        
        phi_12_samples[ii, 1, phi_12_indiv_names] <- c(
          .subset2(phi_12_indiv_prop_list, 1),
          .subset2(phi_12_indiv_prop_list, 2)
        )
      } else {
        if (jj == 1) {
          psi_1_indices[ii, 1, ] <- psi_1_indices[ii - 1, 1, ]
        }
        phi_12_samples[ii, 1, phi_12_indiv_names] <- phi_12_samples[ii - 1, 1, phi_12_indiv_names]
      }
    }
    
    if (ii %% 500 == 0) {
      flog.info(
        sprintf("surv-fit-point-est-3-meld-12--chain-%d: iteration-%d", chain_id, ii), 
        name = base_filename
      )
    }
  }
  
  inner_res <- list(
    psi_1_indices = psi_1_indices,
    phi_12_samples = phi_12_samples,
    psi_2_samples = psi_2_samples
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

flog.info("surv-fit-point-est-3-meld-12: writing to disk", name = base_filename)

saveRDS(
  object = phi_12_samples,
  file = "rds/surv-example/point-est-3-meld-12-phi-12-samples.rds"
)

saveRDS(
  object = psi_2_samples,
  file = "rds/surv-example/point-est-3-meld-12-psi-2-samples.rds"
)

saveRDS(
  object = psi_1_indices,
  file = "rds/surv-example/point-est-3-meld-12-psi-1-indices.rds"
)

flog.info("surv-fit-point-est-3-meld-12: finished writing", name = base_filename)
