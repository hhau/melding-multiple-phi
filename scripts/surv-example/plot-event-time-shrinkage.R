library(tidybayes)
library(magrittr)
library(dplyr)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")
source("scripts/surv-example/GLOBALS.R")

set.seed(sim_seed)

stage_two_event_times <- readRDS("rds/surv-example/stage-two-phi-12-samples.rds")
event_time_names <- dimnames(stage_two_event_times)[[3]]

submodel_one_output <- readRDS("rds/surv-example/submodel-one-output.rds")
stage_one_event_times <- submodel_one_output$samples[, , event_time_names]

stage_one_event_times_tbl <- stage_one_event_times %>% 
  array_to_mcmc_list() %>% 
  gather_draws(event_time[i]) %>% 
  mutate(stage = 1)

stage_two_event_times_tbl <- stage_two_event_times %>% 
  array_to_mcmc_list() %>% 
  gather_draws(event_time[i]) %>% 
  mutate(stage = 2) 

all_tbl <- rbind(stage_one_event_times_tbl, stage_two_event_times_tbl) %>% 
  mutate(
    stage = as.factor(stage), 
    facet_label = factor(
      x = as.character(i),
      levels = i,
      labels = paste0("italic(i) ~ '=' ~", i)
    )
  ) 

# filter down to whom actually had the event:
sim_settings <- readRDS(
  "rds/surv-example/simulation-settings-and-joint-data.rds"
)

index_vec <- which(sim_settings$event_indicator == 1)
event_only_patients <- all_tbl %>% 
  filter(i %in% index_vec) 

# remove all the sub_sampling stuff later if we get the N down
sub_sample <-  c(4, 11, 20, 23, 25, 33)

p_1 <- ggplot(
    data = all_tbl %>% filter(i %in% sub_sample), 
    mapping = aes(x = .value, colour = stage, lty = stage)
  ) +
  geom_density(adjust = 3) +
  facet_wrap(vars(facet_label), scales = "free", labeller = label_parsed) + 
  scale_x_continuous(limits = c(0, 1)) + 
  xlab(expression(italic(t))) + 
  ylab(expression("p"(italic(t)['i'] ~ '|' ~ 'Y'[2] * "," ~ 'Y'[3]))) + # not right - has both stages
  scale_color_manual(values = c('1' = blues[2], '2' = highlight_col)) + 
  labs(colour = "Stage", lty = "Stage")

ggsave_halfheight(
  filename = "plots/surv-example/phi-12-inter-stage-comparison.pdf",
  plot = p_1
)
