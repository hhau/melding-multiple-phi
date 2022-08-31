library(bayesplot)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)
library(stringr)
library(patchwork)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

phi_12_samples <- readRDS("rds/surv-example/stage-two-phi-12-samples.rds")
phi_23_samples <- readRDS("rds/surv-example/stage-two-phi-23-samples.rds")
psi_2_samples <- readRDS("rds/surv-example/stage-two-psi-2-samples.rds")

p_1 <- plot_worst_pars(
  phi_12_samples, 
  c(
    "event_time" = "italic(t)",
    "event_indicator" = "delta"
  ),
)

p_2 <- plot_worst_pars(
  phi_23_samples, 
  facet_name_value_pairs = function(x) {
    res <- str_extract_all(x, "(\\d+,\\d)") %>%
      str_split(",", simplify = TRUE) %>%
      as.numeric() 
    
    sprintf("eta[%d * ',' ~ %d]", res[1], res[2])
  }
)

p_3 <- plot_worst_pars(
  psi_2_samples,
  c(
    "theta_zero" = "theta[0]",
    'theta_one' = "theta[1]",
    "hazard_gamma" = "gamma",
    "alpha" = "alpha"
  )
)

ggsave_halfheight(
  filename = "plots/surv-example/stage-two-phi-12-diags.png",
  plot = p_1
)

ggsave_halfheight(
  filename = "plots/surv-example/stage-two-phi-23-diags.png",
  plot = p_2
)

ggsave_halfheight(
  filename = "plots/surv-example/stage-two-psi-2-diags.png",
  plot = p_3
)
