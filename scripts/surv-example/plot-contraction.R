library(tidyverse)
library(tidybayes)
library(parallel)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

sim_settings <- readRDS(
  file = "rds/surv-example/simulation-settings-and-joint-data.rds"
)

phi_12_stage_one <- readRDS("rds/surv-example/stage-one-phi-12-samples.rds")
phi_12_stage_two <- readRDS("rds/surv-example/stage-two-phi-12-samples.rds")
phi_23_stage_one <- readRDS("rds/surv-example/stage-one-phi-23-samples.rds")
phi_23_stage_two <- readRDS("rds/surv-example/stage-two-phi-23-samples.rds")

phi_12_names <- names(phi_12_stage_one[1, 1, ]) %>%
  grep("event", ., value = TRUE)

phi_23_names <- names(phi_23_stage_one[1, 1, ]) %>%
  grep("^beta", ., value = TRUE)

res <- mclapply(1 : 4, mc.cores = 4, function(job_id) {
  if (job_id == 1) {
    phi_12_stage_one[, , phi_12_names] %>%
      array_to_mcmc_list() %>%
      gather_draws(event_time[i], event_indicator[i]) %>%
      mutate(stage = "1", phi = "12")  
  } else if (job_id == 2) {
    phi_12_stage_two[, , phi_12_names] %>%
      array_to_mcmc_list() %>%
      gather_draws(event_time[i], event_indicator[i]) %>%
      mutate(stage = "2", phi = "12")  
  } else if (job_id == 3) {
    phi_23_stage_one[, , phi_23_names] %>%
      array_to_mcmc_list() %>%
      gather_draws(beta_zero[i]) %>%
      mutate(stage = "1", phi = "23")  
  } else if (job_id == 4) {
    phi_23_stage_two[, , phi_23_names] %>%
      array_to_mcmc_list() %>%
      gather_draws(beta_zero[i]) %>%
      mutate(stage = "2", phi = "23")
  }
}) %>%
  bind_rows()

sub_df <- tibble(
  id = 1 : sim_settings$n_patients,
  indicator = sim_settings$event_indicator
)

plot_df <- res %>%
  select(-c(.chain, .iteration, .draw)) %>%
  group_by(i, .variable, stage) %>%
  summarise(post_sd = sd(.value)) %>%
  pivot_wider(
    names_from = stage, 
    values_from = post_sd, 
    names_prefix = "stage_"
  ) %>%
  left_join(sub_df, by = c("i" = "id")) %>%
  mutate(
    var = factor(
      x = .variable,
      levels = c("event_indicator", "event_time", "beta_zero"),
      labels = c("delta", "italic(t)", "eta")
    ), 
    indicator = as.logical(indicator)
  ) 
  
p1 <- ggplot(plot_df, aes(x = stage_1, y = stage_2, col = indicator)) + 
  geom_point() +
  facet_wrap(vars(var), labeller = label_parsed, scales = "free") +
  geom_abline(slope = 1, intercept = 0) +
  # coord_fixed() +
  xlab("Stage 1 posterior SD") +
  ylab("Stage 2 posterior SD") +
  labs(col = "Had event") +
  guides(col = guide_legend(reverse = TRUE))

ggsave_fullpage(
  filename = "plots/surv-example/phi-inter-stage-posterior-sd.pdf",
  plot = p1,
  adjust_height = -15
)
