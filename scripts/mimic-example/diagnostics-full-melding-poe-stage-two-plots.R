library(stringr)
library(magrittr)

source("scripts/common/plot-settings.R")
source('scripts/common/mcmc-util.R')
source("scripts/common/setup-argparse.R")

parser$add_argument("--stage-two-phi-12-samples")
parser$add_argument("--stage-two-phi-23-samples")
parser$add_argument("--stage-two-psi-2-samples")
parser$add_argument("--phi-12-diagnostic-plot")
parser$add_argument("--phi-23-diagnostic-plot")
parser$add_argument("--output-full-phi-23-trace-plot")
args <- parser$parse_args()

phi_12_samples <- readRDS(args$stage_two_phi_12_samples)
phi_23_samples <- readRDS(args$stage_two_phi_23_samples)
psi_2_samples <- readRDS(args$stage_two_psi_2_samples)

event_time_names <- phi_12_samples[1, 1, ] %>%
  names() %>%
  grep(pattern = 'event_time', x = ., value = T)

n_icustays <- length(event_time_names)

n_segement_per_icustay <- names(phi_23_samples[1, 1, ]) %>%
  grep(pattern = 'eta_slope', x = ., value = T) %>%
  length() %>%
  divide_by(n_icustays)

eta_grid <- expand.grid(icustay_id = 1 : n_icustays, eta_alphabetic = c('b', 'a'))
eta_grid$eta_numeric <- rep(1 : 2, each = n_icustays)

n_theta <- names(psi_2_samples[1, 1, ]) %>%
  grep(pattern = 'theta', x = .) %>%
  length()

p1 <- plot_worst_pars(
  phi_12_samples[, , event_time_names],
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('event_time[%d]', 1 : n_icustays)
    )

    values <- c(
      sprintf("italic(T)[%d]", 1 : n_icustays)
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = args$phi_12_diagnostic_plot,
  plot = p1
)

p1_2 <- mcmc_trace(x = phi_23_samples)

ggsave_base(
  filename = args$output_full_phi_23_trace_plot,
  plot = p1_2,
  height = 48,
  width = 48
)

p2 <- plot_worst_pars(
  phi_23_samples,
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('eta_slope[%d,%d]', eta_grid$icustay_id, eta_grid$eta_numeric),
      sprintf('breakpoint[%d]', 1 : n_icustays)
    )

    values <- c(
      sprintf("eta[1 * ',' ~ %d]^{%s}", eta_grid$icustay_id, eta_grid$eta_alphabetic),
      sprintf('kappa[%d]', 1 : n_icustays)
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = args$phi_23_diagnostic_plot,
  plot = p2
)

p3 <- plot_worst_pars(
  psi_2_samples,
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('theta[%d]', 1 : n_theta),
      'hazard_gamma',
      'dd_gamma',
      'alpha'
    )

    values <- c(
      sprintf('theta[%d]', 1 : n_theta),
      'gamma'['rf'],
      'gamma'['dd'],
      'alpha'
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = args$output,
  plot = p3
)

# fix phi_12 meld phi_23
fp12mp23_psi_2 <- readRDS('rds/mimic-example/stage-two-poe-psi-2-samples-fix-phi-12-meld-phi-23.rds')
fp12mp23_phi_23 <- readRDS('rds/mimic-example/stage-two-poe-phi-23-samples-fix-phi-12-meld-phi-23.rds')

p4 <- plot_worst_pars(
  fp12mp23_psi_2,
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('theta[%d]', 1 : n_theta),
      'hazard_gamma',
      'dd_gamma',
      'alpha'
    )

    values <- c(
      sprintf('theta[%d]', 1 : n_theta),
      'gamma'['rf'],
      'gamma'['dd'],
      'alpha'
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = 'plots/mimic-example/psi-2-stage-2-fix-phi-12-meld-phi-23-diagnostics.png',
  plot = p4
)

p5 <- plot_worst_pars(
  fp12mp23_phi_23,
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('eta_slope[%d,%d]', eta_grid$icustay_id, eta_grid$eta_numeric),
      sprintf('breakpoint[%d]', 1 : n_icustays)
    )

    values <- c(
      sprintf("eta[1 * ',' ~ %d]^{%s}", eta_grid$icustay_id, eta_grid$eta_alphabetic),
      sprintf('kappa[%d]', 1 : n_icustays)
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = 'plots/mimic-example/phi-23-stage-2-fix-phi-12-meld-phi-23-diagnostics.png',
  plot = p5
)

# meld phi_12 fix phi_23
mp12fp23_psi_2 <- readRDS('rds/mimic-example/stage-two-poe-psi-2-samples-meld-phi-12-fix-phi-23.rds')
mp12fp23_phi_12 <- readRDS('rds/mimic-example/stage-two-poe-phi-12-samples-meld-phi-12-fix-phi-23.rds')

p6 <- plot_worst_pars(
  mp12fp23_psi_2,
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('theta[%d]', 1 : n_theta),
      'hazard_gamma',
      'dd_gamma',
      'alpha'
    )

    values <- c(
      sprintf('theta[%d]', 1 : n_theta),
      'gamma'['rf'],
      'gamma'['dd'],
      'alpha'
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = 'plots/mimic-example/psi-2-stage-2-meld-phi-12-fix-phi-23-diagnostics.png',
  plot = p6
)

p7 <- plot_worst_pars(
  mp12fp23_phi_12[, , sprintf('event_time[%d]', 1 : n_icustays)],
  facet_name_value_pairs = function(i) {
    names <- c(
      sprintf('event_time[%d]', 1 : n_icustays)
    )

    values <- c(
      sprintf("italic(T)[%d]", 1 : n_icustays)
    )

    names(values) <- names
    values[i]
  }
)

ggsave_halfheight(
  filename = 'plots/mimic-example/phi-12-stage-2-meld-phi-12-fix-phi-23-diagnostics.png',
  plot = p7
)
