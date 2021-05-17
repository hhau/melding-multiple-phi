source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

population_count_subposterior_samples <- readRDS(
  'rds/owls-example/count-data-subposterior-samples.rds'
)

actual_pars <- c('v[1]', 'v[2]', 'v[6]', 'fec')

p1 <- plot_worst_pars(
 population_count_subposterior_samples[, , actual_pars],
 facet_name_value_pairs = function(x) {
    c(
     'v[1]' = 'alpha[0]',
     'v[2]' = 'alpha[2]', 
     'v[6]' = 'alpha[6]', 
     'fec' = 'rho'
    )[x] 
 }
)

ggsave_halfheight(
  filename = 'plots/owls-example/stage-one-diagnostics-population-count.png',
  plot = p1
)
