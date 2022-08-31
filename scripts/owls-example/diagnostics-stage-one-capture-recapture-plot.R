library(bayesplot)

source("scripts/common/plot-settings.R")

capture_recapture_submodel_samples <- readRDS(
  "rds/owls-example/capture-recapture-subposterior-samples.rds"
)

pars <- sprintf("v[%d]", c(1, 2))

# model 3

stage_one_plot_labeler <- as_labeller(
  x = c(
    'v[1]' = 'alpha[0]',
    'v[2]' = 'alpha[2]'
  ),
  default = label_parsed
)

capture_recapture_traceplot <- mcmc_trace(
  capture_recapture_submodel_samples[, , c(pars), drop = FALSE],
  facet_args = list(labeller = stage_one_plot_labeler)
)

ggsave_halfheight(
  filename = "plots/owls-example/stage-one-diagnostics-capture-recapture.png",
  plot = capture_recapture_traceplot
)

