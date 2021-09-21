source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

library(dplyr)
library(tidybayes)
library(scales)
library(RColorBrewer)
library(magrittr)
library(patchwork)

original_model_samples <- readRDS("rds/owls-example/original-ipm-samples.rds")
capture_recapture_submodel_samples <- readRDS("rds/owls-example/capture-recapture-subposterior-samples.rds")
count_data_submodel_samples <- readRDS("rds/owls-example/count-data-subposterior-samples.rds")
fecunditiy_submodel_samples <- readRDS("rds/owls-example/fecundity-subposterior-samples.rds")
melded_model_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")
normal_approx_samples <- readRDS("rds/owls-example/melded-posterior-normal-approx-samples.rds")
melded_model_log_pooling_phi_samples <- readRDS("rds/owls-example/melded-phi-samples-log-pooling.rds")
melded_model_lin_pooling_phi_samples <- readRDS("rds/owls-example/melded-phi-samples-lin-pooling.rds")

# a tidybayes multiple-interval approach
original_model_tidybayes_tbl <- original_model_samples %>%
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "original-ipm-model")

recode_vec <- c(
  "1" = "alpha[0]",
  "2" = "alpha[2]"
)

original_model_tidybayes_tbl$orig_par <- original_model_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

capture_recapture_tidybayes_tbl <- capture_recapture_submodel_samples %>%
  array_to_mcmc_list() %>%
  gather_draws(v[i]) %>%
  filter(i %in% c(1, 2)) %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "capture-recapture-submodel")

capture_recapture_tidybayes_tbl$orig_par <- capture_recapture_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

count_data_tidybayes_tbl <- count_data_submodel_samples %>%
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "count-data-submodel")

count_data_tidybayes_tbl$orig_par <- count_data_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

melded_model_tidybayes_tbl <- melded_model_samples %>%
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "melded-model-z")  

melded_model_tidybayes_tbl$orig_par <- melded_model_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

normal_approx_tidybayes_tbl <- normal_approx_samples %>%
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "melded-model-a-normal-approx")  

normal_approx_tidybayes_tbl$orig_par <- normal_approx_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

melded_log_pooling_tidybayes_tbl <- melded_model_log_pooling_phi_samples %>% 
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "melded-model-log-pooling")  

melded_log_pooling_tidybayes_tbl$orig_par <- melded_log_pooling_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

melded_lin_pooling_tidybayes_tbl <- melded_model_lin_pooling_phi_samples %>% 
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "melded-model-lin-pooling")  

melded_lin_pooling_tidybayes_tbl$orig_par <- melded_lin_pooling_tidybayes_tbl$i %>%
  recode(!!!recode_vec)

# Try to wrangle the fecundity samples into the same form
# Manually do it
widths_vec <- c(.5, .8, .95, .99)
res <- lapply(widths_vec, function(a_ci_width) {
  all_samples <- as.vector(fecunditiy_submodel_samples)
  tibble(
    i = NA,
    .variable = "fec",
    .value = quantile(all_samples, 0.5),
    .upper = quantile(all_samples, 0.5 + (a_ci_width / 2)),
    .lower = quantile(all_samples, 0.5 - (a_ci_width / 2)),
    .width = a_ci_width,
    .point = "median",
    .interval = "qi",
    model_type = "fecundity-submodel",
    orig_par = "rho"
  )
})

fecundity_tidybayes_tbl <- bind_rows(res)

plot_tbl <- bind_rows(
  original_model_tidybayes_tbl,
  capture_recapture_tidybayes_tbl,
  count_data_tidybayes_tbl,
  fecundity_tidybayes_tbl,
  melded_model_tidybayes_tbl,
  normal_approx_tidybayes_tbl,
  melded_log_pooling_tidybayes_tbl,
  melded_lin_pooling_tidybayes_tbl
)

plot_tbl$orig_par <- plot_tbl$orig_par %>%  
  tidyr::replace_na(
    "rho"
  )

p_2 <- ggplot(
  data = plot_tbl,
  aes(y = model_type, x = .value, colour = model_type)
) +
  facet_wrap(
    vars(orig_par),
    scales = "free",
    ncol = 3,
    nrow = 1,
    labeller = label_parsed
  ) +
  ggdist::geom_interval(
    mapping = aes(xmin = .lower, xmax = .upper),
    alpha = rescale(1 - plot_tbl$.width, to = c(0.1, 1)),
    size = 9,
    orientation = "horizontal"
  ) +
  scale_size_continuous(range = c(12, 20)) +
  scale_color_manual(
    aesthetics = "colour",
    name = "Model",
    values = c(
      "original-ipm-model" = "#000000",
      "capture-recapture-submodel" = highlight_col,
      "count-data-submodel" = greens[2],
      "fecundity-submodel" = "#EE3377",
      "melded-model-z" = blues[1],
      "melded-model-log-pooling" = blues[2],
      "melded-model-lin-pooling" = blues[3],
      "melded-model-a-normal-approx" = "#666666"
    ),
    labels = c(
      "original-ipm-model" = expression("p"["ipm"]),
      "capture-recapture-submodel" = expression("p"[1]),
      "count-data-submodel" = expression("p"[2]),
      "fecundity-submodel" = expression("p"[3]),
      "melded-model-z" = expression("p"["meld," ~ "PoE"]),
      "melded-model-log-pooling" = expression("p"["meld," ~ "log"]),
      "melded-model-lin-pooling" = expression("p"["meld," ~ "lin"]),
      "melded-model-a-normal-approx" = expression(widehat("p")["meld"])
    ),
    guide = guide_legend(
      reverse = TRUE,
      override.aes = list(size = 4)
    )
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = rel(0.9))
  ) +
  xlab("") +
  ylab("") 

ggsave_fullpage(
  filename = "plots/owls-example/subposteriors.pdf",
  plot = p_2,
  adjust_height = -13
)

# this plot filters out the count data submodel and its' wide credible intervals
# all the other models look spot on the same.
p_melding_only <- ggplot(
  data = plot_tbl %>% 
    filter(grepl(pattern = '(original|melded)', x = model_type)),
  aes(y = model_type, x = .value, colour = model_type)
) +
  facet_wrap(
    vars(orig_par),
    scales = "free",
    ncol = 3,
    nrow = 1,
    labeller = label_parsed
  ) +
  ggdist::geom_interval(
    mapping = aes(xmin = .lower, xmax = .upper),
    alpha = plot_tbl %>%
      filter(grepl(pattern = '(original|melded)', x = model_type)) %>%
      mutate(alpha = rescale(1 - .width, to = c(0.2, 1))) %>%
      pull(alpha),
    size = 9,
    orientation = "horizontal"
  ) +
  scale_size_continuous(range = c(12, 20)) +
  scale_color_manual(
    aesthetics = "colour",
    name = "Model",
    values = c(
      "original-ipm-model" = "#000000",
      "capture-recapture-submodel" = highlight_col,
      "count-data-submodel" = greens[2],
      "fecundity-submodel" = "#EE3377",
      "melded-model-z" = blues[1],
      "melded-model-log-pooling" = blues[2],
      "melded-model-lin-pooling" = blues[3],
      "melded-model-a-normal-approx" = "#666666"
    ),
    labels = c(
      "original-ipm-model" = expression("p"["ipm"]),
      "capture-recapture-submodel" = expression("p"[1]),
      "count-data-submodel" = expression("p"[2]),
      "fecundity-submodel" = expression("p"[3]),
      "melded-model-z" = expression("p"["meld," ~ "PoE"]),
      "melded-model-log-pooling" = expression("p"["meld," ~ "log"]),
      "melded-model-lin-pooling" = expression("p"["meld," ~ "lin"]),
      "melded-model-a-normal-approx" = expression(widehat("p")["meld"])
    ),
    guide = guide_legend(
      reverse = TRUE,
      override.aes = list(size = 4)
    )
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = rel(0.9))
  ) +
  xlab("") +
  ylab("")

ggsave_fullpage(
  filename = "plots/owls-example/subposteriors-melding-only.pdf",
  plot = p_melding_only,
  adjust_height = -15
)

p_patchwork <- (p_2 / p_melding_only) + plot_layout(heights = c(7 / 12, 5 / 12))

ggsave_fullpage(
  filename = "plots/owls-example/subposteriors-both-patchwork.pdf",
  plot = p_patchwork,
  adjust_height = -7
)
