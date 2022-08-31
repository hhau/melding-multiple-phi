library(rstan)
library(magrittr)
library(tidybayes)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--fluid-submodel-stan-data")
parser$add_argument("--fluid-piecewise-stan-model")
parser$add_argument("--mimic-globals")
parser$add_argument("--output-array")
parser$add_argument("--output-plot-mu")
args <- parser$parse_args()

source(args$mimic_globals)

stan_data <- readRDS(args$fluid_submodel_stan_data)

init_function <- function(i) {
  n_icu_stays <- stan_data$n_icu_stays

  list(
    y_sigma = abs(rnorm(n = n_icu_stays, mean = 1, sd = 1)) %>%
      array(dim = c(n_icu_stays)),
    eta_zero_raw = rlnorm(
      n = n_icu_stays,
      meanlog = 1.61,
      sdlog = 0.47
    ) %>%
      array(dim = c(n_icu_stays)),
    breakpoint_raw = rbeta(
      n = n_icu_stays,
      shape1 = 5,
      shape2 = 15
    ) %>%
      array(dim = c(n_icu_stays)),
    eta_slope = abs(rnorm(
      n = n_icu_stays * 2,
      mean = 2.5,
      sd = 1.0
    )) %>%
    array(dim = c(n_icu_stays, 2))
  )
}

prefit <- stan_model(args$fluid_piecewise_stan_model)

model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = N_CHAIN,
  chains = N_CHAIN,
  warmup = 1000,
  iter = N_POST_WARMUP_MCMC + 1000,
  init = init_function,
  control = list(
    adapt_delta = 0.98,
    max_treedepth = 12
  )
)

long_samples <- model_fit %>%
  gather_draws(
    y_sigma[i],
    eta_zero[i],
    eta_slope[i, b],
    breakpoint[i],
    breakpoint_raw[i]
  )

array_samples <- model_fit %>%
  as.array(
    pars = c(
      'y_sigma',
      'eta_zero',
      'eta_slope',
      'breakpoint',
      'breakpoint_raw'
    )
  )

plot_mu_samples <- model_fit %>%
  gather_draws(plot_mu[i, p])

saveRDS(file = args$output, object = long_samples)
saveRDS(file = args$output_array, object = array_samples)
saveRDS(file = args$output_plot_mu, object = plot_mu_samples)