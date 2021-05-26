library(rstan)
library(dplyr)
library(parallel)
library(abind)
library(stringr)

source("scripts/common/logger-setup.R")
source("scripts/mimic-example/GLOBALS.R")
source("scripts/common/setup-argparse.R")

parser$add_argument("--pf-event-time-samples-array")
parser$add_argument("--fluid-model-samples-array")
parser$add_argument("--baseline-data")
parser$add_argument("--psi-step-stan-model")
parser$add_argument("--phi-step-indiv-stan-model")
parser$add_argument("--output-phi-12-samples")
parser$add_argument("--output-phi-23-samples")
parser$add_argument("--output-psi-1-indices")
parser$add_argument("--output-psi-3-indices")
args <- parser$parse_args()

submodel_one_output <- readRDS(args$pf_event_time_samples_array)
submodel_three_output <- readRDS(args$fluid_model_samples_array)
submodel_two_data <- readRDS(args$baseline_data)

n_icu_stays <- nrow(submodel_two_data)

stopifnot(
  n_icu_stays == 'eta_zero' %>%
    grep(names(submodel_three_output[1, 1,]), value = TRUE) %>%
    length(),
  n_icu_stays == 'event_time' %>%
    grep(names(submodel_one_output[1, 1, ]), value = TRUE) %>%
    length()
)

submodel_one_samples <- submodel_one_output
n_iter_submodel_one <- dim(submodel_one_samples)[1]
n_chain_submodel_one <- dim(submodel_one_samples)[2]

submodel_three_samples <- submodel_three_output
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
  sprintf("event_time[%d]", 1 : n_icu_stays),
  sprintf("event_indicator[%d]", 1 : n_icu_stays)
)

event_time_names <- grep("event_time", phi_12_names, value = TRUE)
event_indicator_names <- grep("event_indicator", phi_12_names, value = TRUE)

phi_23_names <- c(
  grep("eta_slope", names(submodel_three_samples[1, 1, ]), value = TRUE),
  grep("breakpoint\\[", names(submodel_three_samples[1, 1, ]), value = TRUE)
)

breakpoint_names <- grep('breakpoint', phi_23_names, value = TRUE)
eta_slope_names <- grep('eta_slope', phi_23_names, value = TRUE)

n_phi_12 <- length(phi_12_names)
n_phi_23 <- length(phi_23_names)

n_segments <- grep(
  'eta_slope\\[1,(.+)\\]',
  names(submodel_three_samples[1, 1, ]),
  value = TRUE
) %>%
  length()

center_keep_intercept <- function(x) {
  p <- ncol(x)
  res <- scale(x[, 2 : p], center = TRUE, scale = TRUE)
  cbind(1, res)
}

baseline_data_x_mat <- model.matrix(
  ~ 1 + aniongap_median + bicarbonate_median + creatinine_median + chloride_median +
  glucose_median + hematocrit_median + hemoglobin_median + platelet_median +
  potassium_median + ptt_median + inr_median + pt_median + sodium_median +
  bun_median + wbc_median + gender + age_at_icu_adm,
  data = submodel_two_data %>% as.data.frame()
) %>%
  center_keep_intercept()

n_theta <- ncol(baseline_data_x_mat)

stan_data_psi_base <- list(
  n_icu_stays = n_icu_stays,
  n_theta = n_theta,
  n_segments = n_segments,
  baseline_data_x = baseline_data_x_mat,
  log_crude_event_rate = log(mean(submodel_one_samples[, , event_time_names]))
)

flog.info("MIMIC-stage-two: compiling models", name = base_filename)

prefit_psi_two_step <- stan_model(args$psi_step_stan_model)

# NB: we can use the same stan file/fit for both phi_{1 \cap 2} and
# phi_{2 \cap 3}, because we can control which parameters change / which
# are the same.
# and because we do both phi_{1 \cap 2} and phi_{2 \cap 3} element-at-a-time
prefit_phi_step_indiv <- stan_model(args$phi_step_indiv_stan_model)

stanfit_phi_step_indiv <- sampling(
  prefit_phi_step_indiv,
  data = list(
    n_theta = stan_data_psi_base$n_theta,
    n_segments = stan_data_psi_base$n_segments
  ),
  chains = 1,
  cores = 1,
  iter = 1,
  refresh = 0
)

psi_two_names <- stanfit_phi_step_indiv %>%
  as.array() %>%
  magrittr::extract(1, 1, ) %>%
  names() %>%
  grep(".*(theta|gamma|alpha).*", x = ., value = TRUE)

psi_two_names_simple <- psi_two_names %>%
  str_split('\\[') %>%
  lapply(function(x) x[1]) %>%
  unlist() %>%
  unique()

n_psi_two_pars <- length(psi_two_names)

psi_initialiser <- function(x) {
  list(
    theta = x[sprintf('theta[%d]', 1 : n_theta)],
    hazard_gamma = x['hazard_gamma'],
    dd_gamma = x['dd_gamma'],
    alpha = x['alpha']
  )
}

# general MCMC options -- should be in GLOBALS really
n_stage_two_chain <- N_CHAIN
n_stage_two_iter <- N_POST_WARMUP_MCMC # 4e4 for min tail_ess in phi of ~500

list_res <- mclapply(1 : n_stage_two_chain, mc.cores = N_CHAIN, function(chain_id) {
  # set up containers, remember we are abind'ing over the chains
  # phi - event times (no censoring yet)
  flog.info(
    sprintf("MIMIC-stage-two--chain-%d: initialising", chain_id),
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
    theta = rnorm(n = n_theta, mean = 0.5, sd = 0.2),
    hazard_gamma = abs(rnorm(n = 1, mean = 0.6, sd = 0.2)),
    dd_gamma = rexp(n = 1, rate = 0.75),
    alpha = rnorm(n = 1, mean = -0.2, sd = 0.1)
  )

  flog.info(
    sprintf("MIMIC-stage-two--chain-%d: initialised okay", chain_id),
    name = base_filename
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

    phi_23_curr_list <- list(
      breakpoint = as.numeric(
        phi_23_samples[ii - 1, 1, breakpoint_names]
      ),
      eta_slope = array(
        data = phi_23_samples[ii - 1, 1, eta_slope_names],
        dim = c(n_icu_stays, n_segments)
      )
    )

    # psi step
    psi_step_data <- c(
      stan_data_psi_base,
      phi_12_curr_list,
      phi_23_curr_list
    )

    psi_step <- sampling(
      object = prefit_psi_two_step,
      data = psi_step_data,
      init = list(psi_initialiser(psi_2_samples[ii - 1, 1, ])),
      include = TRUE,
      pars = psi_two_names_simple,
      chains = 1,
      iter = 10,
      warmup = 9,
      refresh = 0
    )

    psi_2_samples[ii, 1, ] <- as.array(psi_step, pars = psi_two_names)
    psi_2_curr_list <- psi_initialiser(psi_2_samples[ii, 1, ])

    # phi_{1 \cap 2} and psi_{1} step -- do 1-at-a-time
    scan_order <- sample.int(n_icu_stays, n_icu_stays, replace = FALSE)
    for (jj in 1 : n_icu_stays) {
      indiv_index <- scan_order[jj]
      phi_12_iter_index_prop <- sample.int(n = n_iter_submodel_one, size = 1)
      phi_12_chain_index_prop <- sample.int(n = n_chain_submodel_one, size = 1)
      phi_12_indiv_names <- phi_12_names[c(indiv_index, indiv_index + n_icu_stays)]

      y2_indiv <- list(
        baseline_data_x = baseline_data_x_mat[indiv_index, ]
      )

      phi_23_indiv_curr_sublist <- list(
        breakpoint = phi_23_curr_list$breakpoint[indiv_index],
        eta_slope = phi_23_curr_list$eta_slope[indiv_index, ]
      )

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
        )
      )

      phi_12_indiv_curr_list <- list(
        event_time = phi_12_samples[ii - 1, 1, phi_12_indiv_names[1]],
        event_indicator = phi_12_samples[ii - 1, 1, phi_12_indiv_names[2]]
      )

      log_prob_phi_12_indiv_prop <- log_prob(
        stanfit_phi_step_indiv,
        upars = unconstrain_pars(
          stanfit_phi_step_indiv,
          pars = c(psi_2_curr_list, y2_indiv, phi_12_indiv_prop_list, phi_23_indiv_curr_sublist)
        )
      )

      log_prob_phi_12_indiv_curr <- log_prob(
        stanfit_phi_step_indiv,
        upars = unconstrain_pars(
          stanfit_phi_step_indiv,
          pars = c(psi_2_curr_list, y2_indiv, phi_12_indiv_curr_list, phi_23_indiv_curr_sublist)
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

    # update phi_12_curr for next step
    phi_12_curr_list <- list(
      event_time = as.numeric(
        phi_12_samples[ii, 1, event_time_names]
      ),
      event_indicator = as.integer(
        phi_12_samples[ii, 1, event_indicator_names]
      )
    )

    # phi_{2 \cap 3} step -- now also do 1 at a time.
    scan_order_2 <- sample.int(n_icu_stays, n_icu_stays, replace = FALSE)
    for (kk in 1 : n_icu_stays) {
      indiv_index <- scan_order_2[kk]
      phi_23_iter_index_prop <- sample.int(n = n_iter_submodel_three, size = 1)
      phi_23_chain_index_prop <- sample.int(n = n_chain_submodel_three, size = 1)
      phi_23_indiv_names <- phi_23_names[
        c(indiv_index, indiv_index + n_icu_stays, indiv_index + (2 * n_icu_stays))
      ]

      y2_indiv <- list(
        baseline_data_x = baseline_data_x_mat[indiv_index, ]
      )

      phi_12_indiv_curr_sublist <- list(
        event_time = phi_12_curr_list$event_time[indiv_index],
        event_indicator = phi_12_curr_list$event_indicator[indiv_index]
      )

      phi_23_indiv_prop_list <- list(
        eta_slope = array(
          data = submodel_three_samples[
            phi_23_iter_index_prop,
            phi_23_chain_index_prop,
            phi_23_indiv_names[1 : 2]
          ],
          dim = c(n_segments)
        ),
        breakpoint = submodel_three_samples[
          phi_23_iter_index_prop,
          phi_23_chain_index_prop,
          phi_23_indiv_names[3]
        ]
      )

      phi_23_indiv_curr_list <- list(
        eta_slope = phi_23_curr_list$eta_slope[indiv_index, ],
        breakpoint = phi_23_curr_list$breakpoint[indiv_index]
      )

      log_prob_phi_23_step_indiv_prop <- log_prob(
        object = stanfit_phi_step_indiv,
        upars = unconstrain_pars(
          stanfit_phi_step_indiv,
          pars = c(
            psi_2_curr_list,
            y2_indiv,
            phi_12_indiv_curr_sublist,
            phi_23_indiv_prop_list
          )
        )
      )

      log_prob_phi_23_step_indiv_curr <- log_prob(
        object = stanfit_phi_step_indiv,
        upars = unconstrain_pars(
          stanfit_phi_step_indiv,
          pars = c(
            psi_2_curr_list,
            y2_indiv,
            phi_12_indiv_curr_sublist,
            phi_23_indiv_curr_list
          )
        )
      )

      log_alpha_23_step_indiv <- log_prob_phi_23_step_indiv_prop - log_prob_phi_23_step_indiv_curr

      if (runif(1) < exp(log_alpha_23_step_indiv)) {
        if (kk == 1) {
          psi_3_indices[ii, 1, ] <- c(phi_23_iter_index_prop, phi_23_chain_index_prop)
        }

        phi_23_samples[ii, 1, phi_23_indiv_names] <- unlist(phi_23_indiv_prop_list)
      } else {
        if (kk == 1) {
          psi_3_indices[ii, 1, ] <- psi_3_indices[ii - 1, 1, ]
        }

        phi_23_samples[ii, 1, phi_23_indiv_names] <-  phi_23_samples[ii - 1, 1, phi_23_indiv_names]
      }

    }

    if (ii %% 500 == 0) {
      flog.info(
        sprintf("MIMIC-stage-two--chain-%d: iteration-%d", chain_id, ii),
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

flog.info("MIMIC-stage-two: writing to disk", name = base_filename)

saveRDS(
  object = phi_12_samples,
  file = args$output_phi_12_samples
)

saveRDS(
  object = phi_23_samples,
  file = args$output_phi_23_samples
)

saveRDS(
  object = psi_2_samples,
  file = args$output
)

saveRDS(
  object = psi_1_indices,
  file = args$output_psi_1_indices
)

saveRDS(
  object = psi_3_indices,
  file = args$output_psi_3_indices
)
