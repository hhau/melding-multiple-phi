library(bayesplot)

source("scripts/common/plot-settings.R")

fecunditiy_submodel_samples <- readRDS(
  "rds/owls-example/fecundity-subposterior-samples.rds"
)

fecundity_traceplot <- mcmc_trace(
  fecunditiy_submodel_samples 
) + 
  ylab(expression(rho)) +
  xlab("Iteration") +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5)
  )

ggsave_halfheight(
  filename = "plots/owls-example/stage-one-diagnostics-fecundity.png",
  plot = fecundity_traceplot
)

