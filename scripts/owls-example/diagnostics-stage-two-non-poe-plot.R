source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

original_model_samples <- readRDS("rds/owls-example/original-ipm-samples.rds")
normal_approx_samples <- readRDS("rds/owls-example/melded-posterior-normal-approx-samples.rds")
melded_model_log_pooling_phi_samples <- readRDS("rds/owls-example/melded-phi-samples-log-pooling.rds")
melded_model_lin_pooling_phi_samples <- readRDS("rds/owls-example/melded-phi-samples-lin-pooling.rds")

pars <- c(sprintf('v[%d]', 1 : 6), 'fec', sprintf('bp[%d]', 1 : 25))
vals <- c(sprintf('alpha[%d]', c(0 : 2, 4 : 6)), 'rho', sprintf('bp[%d]', 1 : 25))
names(vals) <- pars

label_f <- function(x) {
  vals[x]
}

p_orig <- plot_worst_pars(
  original_model_samples[, , pars],
  facet_name_value_pairs = label_f
)

normal_approx_pars <- grep('(fec|im|(^v))', x = pars, value = TRUE)
p_normal <- plot_worst_pars(
  normal_approx_samples[, , normal_approx_pars],
  label_f
)

p_log <- plot_worst_pars(
  melded_model_log_pooling_phi_samples,
  label_f
)

p_lin <- plot_worst_pars(
  melded_model_lin_pooling_phi_samples,
  label_f
)

ggsave_halfheight(
  filename = 'plots/owls-example/diagnostics-original-ipm.png',
  plot = p_orig
)

ggsave_halfheight(
  filename = 'plots/owls-example/stage-two-diagnostics-normal-approx.png',
  plot = p_normal
)

ggsave_halfheight(
  filename = 'plots/owls-example/stage-two-diagnostics-log-pooling.png',
  plot = p_log
)

ggsave_halfheight(
  filename = 'plots/owls-example/stage-two-diagnostics-linear-pooling.png',
  plot = p_lin
)
