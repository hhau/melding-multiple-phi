library(tidybayes)
library(magrittr)
library(dplyr)

source("scripts/common/plot-settings.R")
source("scripts/common/mcmc-util.R")

stage_two_event_times <- readRDS("rds/surv-example/stage-two-phi-samples.rds")
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

all_tbl <- rbind(stage_one_event_times_tbl, stage_two_event_times_tbl)  

# filter down to whom actually had the event:
submodel_one_settings <- readRDS("rds/surv-example/submodel-one-simulation-settings.rds")
index_vec <- which(submodel_one_settings$event_indicator == 1)

all_tbl %>% 
  filter(i %in% index_vec) %>% 
  mutate(
    stage = as.factor(stage),
    i = as.factor(i)
  ) %>% 
  ggplot(aes(x = .value, colour = stage)) +
  geom_density(bw = 0.1) +
  facet_wrap(vars(i), scales = "free") + 
  scale_x_continuous(limits = c(0, 1.2))d
