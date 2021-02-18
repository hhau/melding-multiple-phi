library(tidybayes)
library(dplyr)
library(tidyr)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

stage_one_samples <- readRDS("rds/surv-example/submodel-three-output.rds") %>%
  magrittr::extract2("samples")

stage_two_samples <- readRDS("rds/surv-example/stage-two-phi-23-samples.rds")

param_names <- names(stage_two_samples[1, 1, ])

stage_one_tbl <- stage_one_samples[, , param_names] %>%
  array_to_mcmc_list() %>% 
  gather_draws(beta[i, p]) %>% # will need to change if we change mod 3
  mutate(stage = 1)

stage_two_tbl <- stage_two_samples[, , param_names] %>%
  array_to_mcmc_list() %>%
  gather_draws(beta[i, p]) %>%
  mutate(stage = 2)

plot_tbl <- bind_rows(stage_one_tbl, stage_two_tbl) %>%
  unite("ip", i:p, remove = FALSE) %>%
  mutate(
    stage = as.factor(stage),
    facet_label = factor(
      x = as.character(ip),
      levels = ip,
      labels = paste0("eta[", i, " * ',' ~", p - 1, "]"
      )
    )
  )

sub_sample <- c(4, 11, 20, 23, 25, 33)

p_1 <- ggplot(
  data = plot_tbl %>% filter(i %in% sub_sample),
  mapping = aes(x = .value, colour = stage, lty = stage)
) +
  geom_density() +
  facet_wrap(vars(facet_label), scales = "free", labeller = label_parsed) +
  xlab(expression(eta)) +
  ylab(expression('p'(eta)))  +
  theme(
    axis.text = element_text(size = rel(0.8))
  ) +
  scale_color_manual(values = c('1' = blues[2], '2' = highlight_col)) + 
  labs(colour = "Stage", lty = "Stage")

ggsave_halfheight(
  filename = "plots/surv-example/phi-23-inter-stage-comparison.pdf",
  plot = p_1
)
