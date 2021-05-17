library(rstan)
library(tidybayes)
library(magrittr)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--pf-submodel-stan-data")
parser$add_argument("--pf-submodel-stan-model")
parser$add_argument("--mimic-globals")
parser$add_argument("--output-array")
parser$add_argument("--output-plot-mu")
args <- parser$parse_args()

source(args$mimic_globals)
stan_data <- readRDS(args$pf_submodel_stan_data)
prefit <- stan_model(args$pf_submodel_stan_model)

# relatively expensive, so avoid running an excessive number of iterations
model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = N_CHAIN,
  chains = N_CHAIN,
  warmup = 1000,
  iter = N_POST_WARMUP_MCMC + 1000,
  refresh = 500,
)

long_samples <- model_fit %>%
  gather_draws(y_sigma, beta_zero[i], spline_coef[i, b])

array_samples <- model_fit %>%
  as.array(pars = c('y_sigma', 'beta_zero', 'spline_coef'))

plot_mu_samples <- model_fit %>%
  spread_draws(plot_mu[i, p])

saveRDS(file = args$output, object = long_samples)
saveRDS(file = args$output_array, object = array_samples)
saveRDS(file = args$output_plot_mu, object = plot_mu_samples)
