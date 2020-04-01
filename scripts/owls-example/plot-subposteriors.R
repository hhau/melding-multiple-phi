source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

library(dplyr)
library(tidybayes)
library(scales)
library(RColorBrewer)

# Load:
# - Original model posterior samples (DONE)
# - Capture-Recapture subposterior samples (DONE)
# - Count data subposterior samples (DONE)
# - Stage one samples (TODO)
#   - Maybe don't have to do this if we use the subposterior as the stage-one
#   - Target. 
# - Melded posterior samples (TODO)
original_model_samples <- readRDS("rds/owls-example/original-ipm-samples.rds")
capture_recapture_submodel_samples <- readRDS("rds/owls-example/capture-recapture-subposterior-samples.rds")
count_data_submodel_samples <- readRDS("rds/owls-example/count-data-subposterior-samples.rds")

# try to include rho as well
fecunditiy_submodel_samples <- readRDS("rds/owls-example/fecundity-subposterior-samples.rds")

melded_model_samples <- readRDS("rds/owls-example/melded-posterior-samples.rds")

# a tidybayes multiple-interval approach
original_model_tidybayes_tbl <- original_model_samples %>%
  array_to_mcmc_list() %>%
  gather_draws(v[i], fec) %>%
  filter(i %in% c(1, 2) || .variable == "fec") %>%
  median_qi(.width = c(.5, .8, .95, .99)) %>%
  mutate(model_type = "original-ipm-model")

recode_vec <- c(
  "1" = "alpha[1]",
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
  mutate(model_type = "melded-model")  

melded_model_tidybayes_tbl$orig_par <- melded_model_tidybayes_tbl$i %>%
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
  melded_model_tidybayes_tbl
)

plot_tbl$orig_par <- plot_tbl$orig_par %>%  
  tidyr::replace_na(
    "rho"
  )

p_2 <- ggplot(
  data = plot_tbl,
  aes(y = interaction(model_type, orig_par), x = .value, colour = model_type)
) +
  facet_wrap(
    vars(orig_par),
    scales = "free",
    ncol = 2,
    nrow = 2,
    labeller = label_parsed
  ) +
  geom_intervalh(
    alpha = rescale(1 - plot_tbl$.width, to = c(0.1, 1)),
    size = 9
  ) +
  scale_size_continuous(range = c(12, 20)) +
  scale_color_manual(
    aesthetics = "colour",
    name = "Model",
    values = c(
      "original-ipm-model" = blues[2],
      "capture-recapture-submodel" = highlight_col,
      "count-data-submodel" = greens[2],
      "fecundity-submodel" = "#EE3377",
      "melded-model" = "#666666"
    ),
    labels = c(
      "original-ipm-model" = expression("p"["ipm"]),
      "capture-recapture-submodel" = expression("p"[1]),
      "count-data-submodel" = expression("p"[2]),
      "fecundity-submodel" = expression("p"[3]),
      "melded-model" = expression("p"["meld"])
    ),
    guide = guide_legend(
      reverse = TRUE,
      override.aes = list(size = 4)
    )
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 14)
  ) +
  xlab("") +
  ylab("") 


ggsave_halfheight(
  filename = "plots/owls-example/subposteriors.pdf",
  plot = p_2
)
