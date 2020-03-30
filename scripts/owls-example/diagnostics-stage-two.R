source("scripts/common/plot-settings.R")

library(dplyr)
library(bayesplot)
library(rstan)
library(ggplot2)
library(patchwork)
library(kableExtra)

stage_two_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")
stage_two_diagnostics <- monitor(stage_two_samples)

## numerical diagnostics 
vars <- c("fec", "v[1]", "v[2]")

bulk_min_index <- which.min(stage_two_diagnostics$Bulk_ESS)
tail_min_index <- which.min(stage_two_diagnostics$Tail_ESS)
rhat_max_index <- which.max(stage_two_diagnostics$Rhat)

qois <- names(stage_two_samples[
  1,
  1,
  c(bulk_min_index, tail_min_index, rhat_max_index)
])

diag_vars <- unique(c(vars, qois))

# table output
stage_two_table <- as.data.frame(stage_two_diagnostics)[diag_vars, ] %>%
  tibble::rownames_to_column() %>%
  select(par = rowname, n_eff, Rhat, Bulk_ESS, Tail_ESS) %>%
  kable(
    format = "latex",
    booktabs = TRUE
  ) %>%  
  kable_styling(latex_options = "striped") %>% 
  column_spec(1, "1cm")

cat(
  stage_two_table,
  file = "tex-input/owls-example/appendix-info/0020-stage-two-diagnostics.tex" 
)

## visual diagnostics
p1 <- mcmc_trace(
  x = stage_two_samples[, , vars]
) +
  facet_wrap("parameter", ncol = 1, scales = "free_y") +
  theme(
    legend.position = 'none'
  ) +
  bayesplot:::force_x_axis_in_facets()

p1

p2 <- mcmc_rank_overlay(
  x = stage_two_samples[, , vars],
  ref_line = TRUE
) +
  facet_wrap("parameter", ncol = 1)

p2

res <- p1 + p2 + plot_layout(guides = 'collect')
res

ggsave_fullpage(
  filename = "plots/owls-example/stage-two-diagnostics.png",
  plot = res
)

