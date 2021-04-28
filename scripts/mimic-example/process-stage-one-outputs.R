library(dplyr)

event_time_samples <- readRDS('rds/mimic-example/submodel-1-event-times-samples-long.rds')
fluid_model_samples <- readRDS('rds/mimic-example/submodel-3-fluid-samples-long.rds')

log_crude_event_rate <- event_time_samples %>%
  ungroup() %>% 
  filter(.variable == 'event_time') %>% 
  summarise(log_crude_event_rate = log(mean(.value)))

median_event_time <- event_time_samples %>% 
  group_by(i, .variable) %>% 
  summarise(median = median(.value))

median_fluid_fit_value <- fluid_model_samples %>% 
  filter(.variable %in% c('eta_slope', 'breakpoint')) %>% 
  group_by(.variable, i, b) %>% 
  summarise(
    median = median(.value),
    mean = mean(.value)
  )

saveRDS(
  file = 'rds/mimic-example/median-event-time-data.rds',
  object = list(
    median_event_time,
    log_crude_event_rate
  )
)

saveRDS(
  file = 'rds/mimic-example/median-fluid-fit-data.rds',
  median_fluid_fit_value
)  
