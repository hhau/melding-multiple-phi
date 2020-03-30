source("scripts/common/plot-settings.R")

library(dplyr)
library(ggplot2)
library(tibble)

original_model_samples <- readRDS("rds/owls-example/original-ipm-samples.rds")
stage_two_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")

vars <- c("fec", "v[1]", "v[2]")
plot_quantiles <- seq(from = 0.001, to = 1 - 0.001, by = 0.001)

orig_samples_slim <- array(
  original_model_samples[, , vars],
  dim = c(
    dim(original_model_samples)[1] * dim(original_model_samples)[2],
    3
  ), 
  dimnames = list(
    NULL,
    vars
  )
)

meld_samples_slim <- array(
  stage_two_samples[, , vars],
  dim = c(
    dim(stage_two_samples)[1] * dim(stage_two_samples)[2],
    3
  ), 
  dimnames = list(
    NULL,
    vars
  )
)

res <- lapply(vars, function(a_var) {
  tibble(
    quantile = plot_quantiles,
    orig_quantile = quantile(orig_samples_slim[, a_var], probs = plot_quantiles),
    meld_quantile = quantile(meld_samples_slim[, a_var], probs = plot_quantiles),
    param = a_var
  )
}) %>% 
  bind_rows()

p1 <- ggplot(data = res, aes(x = meld_quantile, y = orig_quantile)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(vars(param), ncol = 1, scales = "free")

ggsave_fullpage(
  filename = "plots/owls-example/orig-meld-qq-compare.pdf",
  plot = p1
)
