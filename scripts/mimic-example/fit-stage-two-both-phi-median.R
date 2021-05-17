library(rstan)
library(dplyr)

source("scripts/common/logger-setup.R")
source("scripts/mimic-example/GLOBALS.R")
source("scripts/common/setup-argparse.R")

flog.info(
  'mimic-example: fitting survival model with fixed (median) values',
  name = base_filename
)

parser$add_argument("--pf-event-time-median")
parser$add_argument("--fluid-model-median")
parser$add_argument("--baseline-data")
parser$add_argument("--psi-step-stan-model")
args <- parser$parse_args()

submodel_one_median_both <- readRDS(args$pf_event_time_median)
submodel_three_median <- readRDS(args$fluid_model_median)
submodel_two_data <- readRDS(args$baseline_data)

submodel_one_median <- submodel_one_median_both[[1]]
log_crude_event_rate <- submodel_one_median_both[[2]]

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

n_icu_stays <- nrow(submodel_two_data)
n_segments <- max(submodel_three_median$b, na.rm = TRUE)

stan_data <- list(
  n_icu_stays = n_icu_stays,
  n_segments = n_segments,
  n_theta = ncol(baseline_data_x_mat),
  baseline_data_x = baseline_data_x_mat,
  log_crude_event_rate = as.numeric(log_crude_event_rate),
  event_indicator = submodel_one_median %>%
    filter(.variable == 'event_indicator') %>%
    pull(median),
  event_time = submodel_one_median %>%
    filter(.variable == 'event_time') %>%
    pull(median),
  breakpoint = submodel_three_median %>%
    filter(.variable == 'breakpoint') %>%
    pull(median),
  eta_slope = submodel_three_median %>%
    filter(.variable == 'eta_slope') %>%
    arrange(b, i) %>%
    pull(median) %>%
    array(dim = c(n_icu_stays, n_segments))
)

prefit <- stan_model(args$psi_step_stan_model)
model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = N_CHAIN,
  chains = N_CHAIN,
  warmup = 1000,
  iter = 1000 + N_POST_WARMUP_MCMC
)

median_samples <- as.array(model_fit)

saveRDS(
  file = args$output,
  object = median_samples
)
