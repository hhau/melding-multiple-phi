library(rstan)
library(dplyr)

source("scripts/common/logger-setup.R")
source("scripts/mimic-example/GLOBALS.R")

submodel_one_median_both <- readRDS("rds/mimic-example/median-event-time-data.rds")
submodel_two_data <- readRDS("rds/mimic-example/baseline-covariate-data.rds")
submodel_three_median <- readRDS("rds/mimic-example/median-fluid-fit-data.rds")

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

prefit <- stan_model("scripts/mimic-example/models/surv-psi-step.stan")
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
  file = 'rds/mimic-example/stage-two-median-inputs-psi-2-samples.rds',
  object = median_samples
)
