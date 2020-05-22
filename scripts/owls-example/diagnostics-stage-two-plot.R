library(bayesplot)
library(rstan)
library(ggplot2)
library(patchwork)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

stage_two_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")
 
vars <- c("fec", "v[1]", "v[2]")

stage_two_plot_labeler <- as_labeller(
  x = c(
    'fec' = 'rho',
    'v[1]' = 'alpha[0]',
    'v[2]' = 'alpha[2]'
  ),
  default = label_parsed
)

# we have to force draw these for patchwork to function correctly
p1 <- mcmc_trace(
  x = stage_two_samples[, , vars]
) +
  facet_wrap(
    "parameter",
    ncol = 1,
    scales = "free_y",
    labeller = stage_two_plot_labeler
  ) +
  xlab("Iteration") + 
  theme(
    legend.position = 'none'
  ) +
  bayesplot:::force_x_axis_in_facets()

p1

p2 <- mcmc_rank_overlay(
  x = stage_two_samples[, , vars],
  ref_line = TRUE,
  n_bins = 20
) +
  facet_wrap("parameter", ncol = 1, labeller = stage_two_plot_labeler)

p2

res <- p1 + p2 + plot_layout(guides = 'collect')
res

ggsave_fullpage(
  filename = "plots/owls-example/stage-two-diagnostics.png",
  plot = res
)

# sometimes the discrete parameters have a tendency to play up
# discrete_pars <- c(
#   sprintf("N1[%d]", 1 : 26),
#   sprintf("Nadimm[%d]", 1 : 26),
#   sprintf("NadSurv[%d]", 1 : 26) 
# )

# discrete_trace <- mcmc_trace(
#   stage_two_samples[, , discrete_pars]
# )

# ggsave(
#   filename = "plots/owls-example/stage-two-diagnostics-discrete.png",
#   plot = discrete_trace,
#   height = 30,
#   width = 30
# )
