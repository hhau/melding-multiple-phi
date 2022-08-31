library(stringr)
library(magrittr)

source('scripts/common/plot-settings.R')
source('scripts/common/mcmc-util.R')
source("scripts/common/setup-argparse.R")

parser$add_argument("--both-median-psi-2-samples")
args <- parser$parse_args()

psi_2_samples <- readRDS(args$both_median_psi_2_samples)

n_theta <- names(psi_2_samples[1, 1, ]) %>%
  grep(pattern = 'theta', x = .) %>%
  length()

samples_no_lp__ <-  names(psi_2_samples[1, 1, ]) %>%
  str_detect(pattern = 'lp__', negate = TRUE) %>%
  extract(psi_2_samples, , , .)

p1 <- plot_worst_pars(
  samples_no_lp__,
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
  plot = p1
)
