library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(pbapply)
library(latex2exp)
library(matrixStats)

source("scripts/common/plot-settings.R")
source("scripts/common/logger-setup.R")
source('scripts/common/setup-argparse.R')

parser$add_argument("--baseline-covariate-data")
parser$add_argument("--stage-two-poe-psi-2-samples")
parser$add_argument("--stage-two-poe-phi-23-samples")
parser$add_argument("--stage-two-poe-phi-12-samples")
parser$add_argument("--stage-two-median-inputs-psi-2-samples")
parser$add_argument("--stage-one-phi-12-point-est")
parser$add_argument("--stage-one-phi-23-point-est")
args <- parser$parse_args()

source("scripts/mimic-example/GLOBALS.R")

submodel_two_data <- readRDS(args$baseline_covariate_data)
psi_2_poe_samples <- readRDS(args$stage_two_poe_psi_2_samples)
phi_23_poe_samples <- readRDS(args$stage_two_poe_phi_23_samples)
phi_12_poe_samples <- readRDS(args$stage_two_poe_phi_12_samples)
psi_2_fixed_samples <- readRDS(args$stage_two_median_inputs_psi_2_samples)
phi_12_point_est <- readRDS(args$stage_one_phi_12_point_est)
phi_23_point_est <- readRDS(args$stage_one_phi_23_point_est)

psi_2_names <- dimnames(psi_2_poe_samples)[[3]]

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
  center_keep_intercept() %>%
  as_tibble()

breakpoint_names <- grep("breakpoint", names(phi_23_poe_samples[1, 1,]), value = TRUE)
eta_one_before_names <- grep("*,1]", names(phi_23_poe_samples[1, 1,]), value = TRUE)
eta_one_after_names <- grep("*,2]", names(phi_23_poe_samples[1, 1,]), value = TRUE)

n_plot <- 250
plot_t <- seq(from = 0, to = max(phi_12_poe_samples), length.out = n_plot)
n_iter <- dim(phi_12_poe_samples)[1]
n_thin_per_chain <- 100
n_chain <- 5
thin_vec <- seq.int(from = 1 , to = n_iter, length.out = n_thin_per_chain) %>%
  round()

phi_12_thin_vec <- seq.int(1, n_iter, length.out = 5 * n_thin_per_chain) %>%
  round()

event_time_names <- names(phi_12_poe_samples[1, 1,]) %>%
  grep("time", x = ., value = T)

event_indicator_names <- names(phi_12_poe_samples[1, 1, ]) %>%
  grep("indicator", x = ., value = T)

n_icustays <- length(event_time_names)

# it would be nice to speed this up some, or at least cache the results
res <- pblapply(phi_12_thin_vec, cl = 6, function(iter_id) {
  df <- data.frame(
    event_time = phi_12_poe_samples[iter_id, 1, event_time_names],
    event_indicator = phi_12_poe_samples[iter_id, 1, event_indicator_names],
    iter_id = iter_id
  )
  surv_obj <- survfit(Surv(event_time, event_indicator) ~ 1, data = df) %>%
    surv_summary() %>%
    as_tibble() %>%
    select(-c(std.err, upper, lower)) %>%
    mutate(iter_id = iter_id)

  base_df <- surv_obj[1, ] %>%
    mutate(time = 0, surv = 1, n.event = 0)

  res <- rbind(base_df, surv_obj)
})

df_fixed <- data.frame(
  event_time = phi_12_point_est[[1]] %>%
    filter(.variable == 'event_time') %>%
    pull(median),
  event_indicator = phi_12_point_est[[1]] %>%
    filter(.variable == 'event_indicator') %>%
    pull(median),
  iter_id = max(phi_12_thin_vec) + 1
)

surv_obj_fixed <- survfit(Surv(event_time, event_indicator) ~ 1, data = df_fixed) %>%
  surv_summary() %>%
  as_tibble() %>%
  select(-c(std.err, upper, lower)) %>%
  mutate(iter_id = max(phi_12_thin_vec) + 1)

base_df_fixed <- surv_obj_fixed[1, ] %>%
  mutate(time = 0, surv = 1, n.event = 0)

res_fixed <- rbind(base_df_fixed, surv_obj_fixed) %>%
  mutate(origin = 'stage_one_median_phi_12', origin_alpha = 0.8, origin_step = 'lty_median')

plot_df <- bind_rows(res) %>%
  mutate(origin = 'stage_two_poe_phi_12', origin_alpha = 0.025, origin_step = 'lty_poe') %>%
  bind_rows(res_fixed)

censor_df <- plot_df %>%
  filter(n.censor > 0)

base_plot <- ggplot(
  data = plot_df,
  aes(x = time, y = surv, group = iter_id, col = origin)
) +
  geom_step(
    aes(alpha = origin_alpha, linetype = origin_step, size = origin_step)
  ) +
  geom_point(
    data = censor_df,
    inherit.aes = FALSE,
    mapping = aes(x = time, y = surv, col = origin, alpha = origin_alpha),
    shape = "|",
    size = 3
  ) +
  scale_alpha(guide = FALSE) +
  scale_linetype_manual(
    guide = FALSE,
    values = c(
      lty_median = '3131',
      lty_poe = 'solid'
    )
  ) +
  scale_size_manual(
    guide = FALSE,
    values = c(
      lty_median = 1,
      lty_poe = 0.125
    )
  )

log_cumulative_hazard_time_before_breakpoint <- function(
  t,
  hazard_gamma,
  baseline_covariates,
  theta,
  alpha,
  eta_before
) {
  log_res <- hazard_gamma * log(t) +
    as.numeric(baseline_covariates %*% theta) +
    alpha * eta_before

  return(log_res)
}

log_cumulative_hazard_time_after_breakpoint <- function(
  t,
  breakpoint,
  hazard_gamma,
  baseline_covariates,
  theta,
  alpha,
  eta_before,
  eta_after
) {
  log_res_base <- as.numeric(baseline_covariates %*% theta)
  log_t1 <- hazard_gamma * log(breakpoint) + (alpha * eta_before)
  log_t2 <- (alpha * eta_after) + log(
    exp(hazard_gamma * log(t)) - exp(hazard_gamma * log(breakpoint))
  )

  log_t3 <- logSumExp(c(log_t1, log_t2))
  log_res <- log_res_base + log_t3

  return(log_res)
}

surv_prob <- function(
  t,
  breakpoint,
  hazard_gamma,
  baseline_covariates,
  theta,
  alpha,
  eta_before,
  eta_after
) {
  if (t < breakpoint) {
    log_cumulative_hazard <- log_cumulative_hazard_time_before_breakpoint(
      t,
      hazard_gamma,
      baseline_covariates,
      theta,
      alpha,
      eta_before
    )
  } else {
    log_cumulative_hazard <- log_cumulative_hazard_time_after_breakpoint(
      t,
      breakpoint,
      hazard_gamma,
      baseline_covariates,
      theta,
      alpha,
      eta_before,
      eta_after
    )
  }

  return(exp(-exp(log_cumulative_hazard)))
}

psi_2_thinned_samples_tbl <- matrix(
  psi_2_poe_samples[thin_vec, , ],
  nrow = n_thin_per_chain * n_chain,
  ncol = dim(psi_2_poe_samples)[3],
  dimnames = list(NULL, names(psi_2_poe_samples[1, 1, ]))
) %>%
  as_tibble() %>%
  mutate(
    chain_id = rep(1 : n_chain, each = n_thin_per_chain),
    sample_id = rep(thin_vec, times = n_chain)
  )

phi_23_thinned_samples_tbl <- tibble(
  patient_id = rep(1 : n_icustays, each = n_thin_per_chain * n_chain),
  chain_id = rep(rep(1 : n_chain, each = n_thin_per_chain), times = n_icustays),
  sample_id = rep(thin_vec, times = n_chain * n_icustays),
  breakpoint = phi_23_poe_samples[thin_vec, , breakpoint_names] %>% as.numeric(),
  eta_before = phi_23_poe_samples[thin_vec, , eta_one_before_names] %>% as.numeric(),
  eta_after = phi_23_poe_samples[thin_vec, , eta_one_after_names] %>% as.numeric()
)

surv_prob_samples <- phi_23_thinned_samples_tbl %>%
  left_join(psi_2_thinned_samples_tbl, by = c('sample_id', 'chain_id')) %>%
  left_join(
    baseline_data_x_mat %>%
      mutate(patient_id = 1 : n()),
    by = 'patient_id'
  )

haz_df <- pblapply(plot_t, cl = 6, function(a_t) {
  surv_prob_samples %>%
    rowwise() %>%
    mutate(
      plot_t = a_t,
      surv_prob = surv_prob(
        t = a_t,
        breakpoint = breakpoint,
        hazard_gamma = hazard_gamma,
        baseline_covariates = c(
          V1, aniongap_median, bicarbonate_median, creatinine_median,
          chloride_median, glucose_median, hematocrit_median, hemoglobin_median,
          platelet_median, potassium_median, ptt_median, inr_median,
          pt_median, sodium_median, bun_median, wbc_median,
          genderM, age_at_icu_adm
        ),
        theta = c(
          `theta[1]`, `theta[2]`, `theta[3]`, `theta[4]`, `theta[5]`,
          `theta[6]`, `theta[7]`, `theta[8]`, `theta[9]`, `theta[10]`,
          `theta[11]`, `theta[12]`, `theta[13]`, `theta[14]`, `theta[15]`,
          `theta[16]`, `theta[17]`, `theta[18]`
        ),
        alpha = alpha,
        eta_before = eta_before,
        eta_after = eta_after
      )
    )
}) %>%
  bind_rows()

# important that we average over individuals within each iteration first.
haz_plot_tbl <- haz_df %>%
  select(patient_id, sample_id, plot_t, surv_prob) %>%
  group_by(sample_id, plot_t) %>%
  summarise(surv_prob = mean(surv_prob)) %>%
  ungroup() %>%
  group_by(plot_t) %>%
  summarise(
    mean = mean(surv_prob),
    median = median(surv_prob),
    lower = quantile(surv_prob, 0.9),
    upper = quantile(surv_prob, 0.1)
  ) %>%
  mutate(method = "chained")

with_melded_interval <- base_plot +
  geom_line(
    data = haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, y = mean),
    col = blues[2]
  ) +
  geom_ribbon(
    data = haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, ymin = lower, ymax = upper),
    alpha = 0.2,
    col = blues[2]
  )

## now add the fully fixed version
psi_2_fixed_thinned_samples_tbl <- matrix(
  psi_2_fixed_samples[thin_vec, , psi_2_names],
  nrow = n_thin_per_chain * n_chain,
  ncol = length(psi_2_names),
  dimnames = list(NULL, psi_2_names)
) %>%
  as_tibble() %>%
  mutate(
    chain_id = rep(1 : n_chain, each = n_thin_per_chain),
    sample_id = rep(thin_vec, times = n_chain)
  )

phi_23_fixed_thinned_samples_tbl <- tibble(
  patient_id = rep(1 : n_icustays, each = n_thin_per_chain * n_chain),
  chain_id = rep(rep(1 : n_chain, each = n_thin_per_chain), times = n_icustays),
  sample_id = rep(thin_vec, times = n_chain * n_icustays),
  breakpoint = phi_23_point_est %>%
    filter(.variable == 'breakpoint') %>%
    pull(mean) %>%
    rep(each = n_thin_per_chain * n_chain),
  eta_before = phi_23_point_est %>%
    filter(.variable == 'eta_slope', b == 1) %>%
    pull(mean) %>%
    rep(each = n_thin_per_chain * n_chain),
  eta_after = phi_23_point_est %>%
    filter(.variable == 'eta_slope', b == 2) %>%
    pull(mean) %>%
    rep(each = n_thin_per_chain * n_chain)
)

surv_prob_fixed_samples <- phi_23_fixed_thinned_samples_tbl %>%
  left_join(
    psi_2_fixed_thinned_samples_tbl,
    by = c("chain_id", "sample_id")
  ) %>%
  left_join(
    baseline_data_x_mat %>%
      mutate(patient_id = 1 : n()),
    by = 'patient_id'
  )

haz_fixed_df <- pblapply(plot_t, cl = 6, function(a_t) {
  surv_prob_fixed_samples %>%
    rowwise() %>%
    mutate(
      plot_t = a_t,
      surv_prob = surv_prob(
        t = a_t,
        breakpoint = breakpoint,
        hazard_gamma = hazard_gamma,
        baseline_covariates = c(
          V1, aniongap_median, bicarbonate_median, creatinine_median,
          chloride_median, glucose_median, hematocrit_median, hemoglobin_median,
          platelet_median, potassium_median, ptt_median, inr_median,
          pt_median, sodium_median, bun_median, wbc_median,
          genderM, age_at_icu_adm
        ),
        theta = c(
          `theta[1]`, `theta[2]`, `theta[3]`, `theta[4]`, `theta[5]`,
          `theta[6]`, `theta[7]`, `theta[8]`, `theta[9]`, `theta[10]`,
          `theta[11]`, `theta[12]`, `theta[13]`, `theta[14]`, `theta[15]`,
          `theta[16]`, `theta[17]`, `theta[18]`
        ),
        alpha = alpha,
        eta_before = eta_before,
        eta_after = eta_after
      )
    )
}) %>%
  bind_rows()

haz_fixed_plot_tbl <- haz_fixed_df %>%
  select(patient_id, sample_id, plot_t, surv_prob) %>%
  group_by(sample_id, plot_t) %>%
  summarise(surv_prob = mean(surv_prob)) %>%
  ungroup() %>%
  group_by(plot_t) %>%
  summarise(
    mean = mean(surv_prob),
    median = median(surv_prob),
    lower = quantile(surv_prob, 0.9),
    upper = quantile(surv_prob, 0.1)
  ) %>%
  mutate(method = "fixed")

both_haz_plot_tbl <- bind_rows(haz_plot_tbl, haz_fixed_plot_tbl)

p1 <- base_plot +
  geom_line(
    data = both_haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, y = mean, col = method)
  ) +
  geom_ribbon(
    data = both_haz_plot_tbl,
    inherit.aes = FALSE,
    aes(x = plot_t, ymin = lower, ymax = upper, col = method),
    alpha = 0.2
  ) +
  scale_colour_manual(
    guide = FALSE,
    values = c(
      stage_one_median_phi_12 = blues[2],
      stage_two_poe_phi_12 = highlight_col,
      chained = highlight_col,
      fixed = blues[2]
    ),
    labels = list(
      stage_one_median_phi_12 = TeX('Stage one subposterior median $\\phi_{1 \\bigcap 2}$'),
      stage_two_poe_phi_12 = TeX('Stage two melded posterior $\\phi_{1 \\bigcap 2}$'),
      chained = "Chained Melding - PoE",
      fixed = TeX("Fix $\\phi_{1 \\bigcap 2}$ and $\\phi_{2 \\bigcap 3}$")
    )
  ) +
  labs(colour = "Posterior type") +
  ylab(expression("S"(italic(t)))) +
  xlab(expression(italic(t) ~ (plain("Days since ICU admission")))) +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave_halfheight(
  filename = args$output,
  plot = p1
)
