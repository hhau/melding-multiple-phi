library(tidyverse)
library(tidybayes)
library(parallel)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")
source("scripts/surv-example/GLOBALS.R")

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
      gather_draws(beta[i, p]) %>%
      mutate(stage = "1", phi = "23")  
  } else if (job_id == 4) {
    phi_23_stage_two[, , phi_23_names] %>%
      array_to_mcmc_list() %>%
      gather_draws(beta[i, p]) %>%
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
  group_by(i, p, .variable, stage) %>%
  summarise(post_sd = sd(.value)) %>%
  pivot_wider(
    names_from = stage, 
    values_from = post_sd, 
    names_prefix = "stage_"
  ) %>%
  left_join(sub_df, by = c("i" = "id")) %>%
  unite("plot_var", .variable:p, na.rm = TRUE) %>%
  mutate(
    var = factor(
      x = plot_var,
      levels = c(
        "event_indicator", 
        "event_time", 
        sprintf("beta_%d", 1 : n_long_beta)
      ),
      labels = c(
        "delta[i]", 
        "t[i]", 
        sprintf("eta[%d * ',' ~ 'i']", (1 : n_long_beta) - 1)
      )
    ), 
    indicator = as.logical(indicator)
  ) 
  
p1 <- ggplot(plot_df, aes(x = stage_1, y = stage_2, col = indicator, shape = indicator)) + 
  geom_point() +
  facet_wrap(vars(var), labeller = label_parsed, scales = "free") +
  geom_abline(slope = 1, intercept = 0) +
  # coord_fixed() +
  xlab("Stage 1 posterior SD") +
  ylab("Stage 2 posterior SD") +
  labs(col = "Had event", shape = "Had event") +
  guides(
    col = guide_legend(reverse = TRUE), 
    shape = guide_legend(reverse = TRUE)
  ) +
  scale_colour_manual(
    values = c("TRUE" = blues[2], "FALSE" = highlight_col)
  ) + 
  scale_shape_manual(
    values = c("TRUE" = 3, "FALSE" = 8)
  )

ggsave_fullpage(
  filename = "plots/surv-example/phi-inter-stage-posterior-sd.pdf",
  plot = p1,
  adjust_height = -15
)
