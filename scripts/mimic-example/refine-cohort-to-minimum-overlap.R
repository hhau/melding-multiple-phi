library(dplyr)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--pf-cohort-and-data")
parser$add_argument("--raw-fluid-data")
args <- parser$parse_args()

flog.info(
  'mimic-example: reading data',
  name = base_filename
)

pf_data <- readRDS(args$pf_cohort_and_data)
fluid_data <- readRDS(args$raw_fluid_data)

pf_time_extrema <- pf_data %>%
  group_by(icustay_id) %>%
  summarise(
    pf_min_time = min(time_since_icu_adm),
    pf_max_time = max(time_since_icu_adm)
  )

fluid_time_extrema <- fluid_data %>%
  filter(icustay_id %in% unique(pf_data$icustay_id)) %>%
  group_by(icustay_id) %>%
  summarise(
    fluid_min_time = min(time_since_icu_adm),
    fluid_max_time = max(time_since_icu_adm)
  )

all_data <- inner_join(
  pf_time_extrema,
  fluid_time_extrema,
  by = 'icustay_id'
) %>%
  arrange(icustay_id) %>%
  filter(complete.cases(.)) %>%
  rowwise() %>%
  mutate(
    overlap = max(
      0,
      min(pf_max_time, fluid_max_time) - max(pf_min_time, fluid_min_time)
    ),
    overlap_as_proportion_of_possible = overlap / (
      max(pf_max_time, fluid_max_time) - min(pf_min_time, fluid_min_time)
    )
  )

overlap_proportion <- 0.9

flog.info(
  sprintf('mimic-example: filtering to overlap of %f', overlap_proportion),
  name = base_filename
)

overlap_cohort <- all_data %>%
  filter(overlap_as_proportion_of_possible > overlap_proportion)

plot_tbl <- pf_data %>%
  filter(icustay_id %in% overlap_cohort$icustay_id) %>%
  select(
    -c(subject_id, hadm_id, charttime, specimen_pred, intime),
    value = pf
  ) %>%
  mutate(value_type = 'pf')

both_plot_tbl <- fluid_data %>%
  filter(icustay_id %in% overlap_cohort$icustay_id) %>%
  select(-c(origin, hadm_id)) %>%
  mutate(value_type = 'fluids') %>%
  bind_rows(plot_tbl) %>%
  select(-c(charttime, intime)) %>%
  arrange(icustay_id, time_since_icu_adm)

saveRDS(
  file = args$output,
  object = both_plot_tbl
)
