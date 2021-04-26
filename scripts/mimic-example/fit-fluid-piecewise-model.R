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

init_gen <- function(chain_id) {
  n_icu_stays <- stan_data$n_icu_stays

  list(
    breakpoint_raw = array(
      data = rbeta(n = n_icu_stays, 10, 10),
      dim = n_icu_stays
    ),
    beta_slope = array(
      data = rnorm(n = n_icu_stays * 2, mean = c(5000, 3000), sd = 500)
        %>% abs(),
      dim = c(n_icu_stays, 2)
    ),
    beta_zero = array(
      data = rnorm(n = n_icu_stays, mean = 7500, sd = 500),
      dim = n_icu_stays
    ),
    y_sigma = rnorm(n = 1, mean = 800, sd = 200)
      %>% abs()
  )
}

prefit <- stan_model(args$fluid_piecewise_stan_model)

model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = 5,
  chains = 5,
  init = init_gen,
  control = list(
    adapt_delta = 0.9,
    max_treedepth = 11
  )
)

long_samples <- model_fit %>%
  gather_draws(
    y_sigma,
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