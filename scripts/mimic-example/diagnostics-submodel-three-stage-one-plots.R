library(stringr)
library(magrittr)

source("scripts/common/plot-settings.R")
source('scripts/common/mcmc-util.R')
source("scripts/common/setup-argparse.R")

parser$add_argument("--fluid-submodel-samples-array")
args <- parser$parse_args()

submodel_3_stage_1_samples <- readRDS(args$fluid_submodel_samples_array)

n_icustays <- names(submodel_3_stage_1_samples[1, 1, ]) %>%
  grep(pattern = 'eta_zero', x = ., value = T) %>%
  length()

n_segement_per_icustay <- names(submodel_3_stage_1_samples[1, 1, ]) %>%
  grep(pattern = 'eta_slope', x = ., value = T) %>%
  length() %>%
  divide_by(n_icustays)

eta_grid <- expand.grid(icustay_id = 1 : n_icustays, eta_alphabetic = c('b', 'a'))
eta_grid$eta_numeric <- rep(1 : 2, each = n_icustays)

samples_no_raw <- names(submodel_3_stage_1_samples[1, 1, ]) %>%
  str_detect(pattern = 'raw', negate = TRUE) %>%
  extract(submodel_3_stage_1_samples, , , .)

p1 <- plot_worst_pars(
  samples_no_raw,
  facet_name_value_pairs = function(i) {
    names <- c(
      'y_sigma',
      sprintf('eta_zero[%d]', 1 : n_icustays),
      sprintf('eta_slope[%d,%d]', eta_grid$icustay_id, eta_grid$eta_numeric),
      sprintf('breakpoint[%d]', 1 : n_icustays)
    )

    values <- c(
      'sigma[x]',
      sprintf("eta[0 * ',' ~ %d]", 1 : n_icustays),
      sprintf("eta[1 * ',' ~ %d]^{%s}", eta_grid$icustay_id, eta_grid$eta_alphabetic),
      sprintf('kappa[%d]', 1 : n_icustays)
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = args$output,
  plot = p1
)
