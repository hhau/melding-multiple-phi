library(dplyr)

source('scripts/common/setup-argparse.R')
source('scripts/common/logger-setup.R')

parser$add_argument("--pf-event-time-samples-long")
parser$add_argument("--fluid-model-samples-long")
parser$add_argument("--output-fluid")
args <- parser$parse_args()

flog.info(
  "mimic-example: getting stage one subposterior medians",
  name = base_filename
)

event_time_samples <- readRDS(args$pf_event_time_samples_long)
fluid_model_samples <- readRDS(args$fluid_model_samples_long)

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
  file = args$output,
  object = list(
    median_event_time,
    log_crude_event_rate
  )
)

saveRDS(
  file = args$output_fluid,
  median_fluid_fit_value
)
