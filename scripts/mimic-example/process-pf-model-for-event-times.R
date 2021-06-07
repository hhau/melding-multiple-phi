library(dplyr)
library(splines2)
library(rootSolve)
library(abind)
library(parallel)
library(tidybayes)

source('scripts/common/mcmc-util.R')
source("scripts/common/setup-argparse.R")

parser$add_argument("--pf-data-list-format")
parser$add_argument("--pf-submodel-samples-long")
parser$add_argument("--output-long")
args <- parser$parse_args()

list_data <- readRDS(args$pf_data_list_format)
samples_long <- readRDS(args$pf_submodel_samples_long)

n_icu_stays <- length(list_data)
n_chain <- max(samples_long$.chain)

res <- mclapply(1 : n_chain, mc.cores = n_chain, function(chain_id) {
  chain_data <- samples_long %>%
    filter(.chain == chain_id)

  n_iteration <- max(chain_data$.iteration)
  n_spline_coef <- chain_data %>%
    pull(b) %>%
    n_distinct(na.rm = TRUE)

  beta_zero_samples <- chain_data %>%
    filter(.variable == 'beta_zero') %>%
    pull(.value) %>%
    array(dim = c(n_iteration, n_icu_stays))

  spline_coef_samples <- chain_data %>%
    filter(.variable == 'spline_coef') %>%
    arrange(b, i, .iteration) %>%
    pull(.value) %>%
    array(dim = c(n_iteration, n_icu_stays, n_spline_coef))

  event_samples <- array(
    data = NA,
    dim = c(n_iteration, 1, 2 * n_icu_stays),
    dimnames = list(
      sprintf('iteration_%d', 1 : n_iteration),
      sprintf('chain_%d', chain_id),
      c(
        sprintf('event_time[%d]', 1 : n_icu_stays),
        sprintf('event_indicator[%d]', 1 : n_icu_stays)
      )
    )
  )

  for (iteration_index in 1 : n_iteration) {
    for (icustay_index in 1 : n_icu_stays) {
      beta_zero <- beta_zero_samples[iteration_index, icustay_index]
      spline_coef <- spline_coef_samples[iteration_index, icustay_index, ]

      boundary_knots <- list_data[[icustay_index]]$boundary_knots
      internal_knots <- list_data[[icustay_index]]$internal_knots
      scaled_threshold <- list_data[[icustay_index]]$threshold_ctr
      init <- boundary_knots[1] + 1e-4

      eval_spline <- function(x) {
        sp_mat <- bSpline(x, knots = internal_knots, Boundary.knots = boundary_knots)
        res <- beta_zero + sp_mat %*% spline_coef
        return(res)
      }

      root_f <- function(x) {
        eval_spline(x) - scaled_threshold
      }

      # need to do some error handling for no root splines
      res_rootsolve <- rootSolve::uniroot.all(
        f = root_f,
        interval = c(
          max(boundary_knots[1], 0),
          boundary_knots[2]
        )
      )

      if (length(res_rootsolve) == 0) {
        # if there is no root, then we don't observe the respiratory failure
        # event and they are instead assumed to be discharged, or expire, at
        # their final observation.
        event_time <- boundary_knots[2]
        event_indicator <- 2
      } else {
        event_time <- min(res_rootsolve)

        if (event_time == 0) {
          event_time = 1e-8 # otherwise the beta distribution falls over
        }

        # this is an annoying numeric edge case
        if (event_time == boundary_knots[2]) {
          event_indicator <- 2
        } else {
          event_indicator <- 1
        }
      }

      event_samples[
        iteration_index,
        1,
        c(icustay_index, n_icu_stays + icustay_index)
      ] <- c(event_time, event_indicator)
    }
  }

  return(event_samples)

}) %>%
  abind(along = 2)

saveRDS(
  file = args$output,
  object = res
)

res_long <- res %>%
  array_to_mcmc_list() %>%
  gather_draws(event_time[i], event_indicator[i])

saveRDS(
  file = args$output_long,
  object = res_long
)
