library(stringr)
library(magrittr)

source("scripts/common/plot-settings.R")
source('scripts/common/mcmc-util.R')
source("scripts/common/setup-argparse.R")

parser$add_argument("--pf-submodel-samples-array")
args <- parser$parse_args()

submodel_1_stage_1_samples <- readRDS(args$pf_submodel_samples_array)

n_icustays <- names(submodel_1_stage_1_samples[1, 1, ]) %>%
  grep(pattern = 'beta_zero', x = ., value = T) %>%
  length()

n_spline_coef_per_icustay <- names(submodel_1_stage_1_samples[1, 1, ]) %>%
  grep(pattern = 'spline_coef', x = ., value = T) %>%
  length() %>%
  divide_by(n_icustays)

spline_grid <- expand.grid(1 : n_icustays, 1 : n_spline_coef_per_icustay)

p1 <- plot_worst_pars(
  submodel_1_stage_1_samples,
  facet_name_value_pairs = function(i) {
    names <- c(
      'y_sigma',
      sprintf('beta_zero[%d]', 1 : n_icustays),
      sprintf('spline_coef[%d,%d]', spline_grid$Var1, spline_grid$Var2)
    )

    values <- c(
      'omega',
      sprintf("beta[0 * ',' ~ %d]", 1 : n_icustays),
      sprintf("zeta[%d * ',' ~ %d]", spline_grid$Var1, spline_grid$Var2)
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = args$output,
  plot = p1
)
