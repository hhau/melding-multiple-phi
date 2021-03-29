library(RPostgres)
library(dplyr)
library(readr)
library(withr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(mgcv)

mimic <- dbConnect(
  RPostgres::Postgres(),
  dbname = "mimic",
  host = "localhost",
  port = 5432,
  user = "amanderson",
  password = "postgres",
  timezone = "US/Eastern"
)

icustays <- dbGetQuery(
  mimic,
  statement = 
    "SELECT 
    subject_id,
    hadm_id,
    icustay_id,
    intime
  FROM 
    mimiciii.icustays
  "
) %>% 
  as_tibble()

# these are all inputs, using the standard mimic queries
inputs_cv <- dbGetQuery(
  mimic,
  statement = read_file('scripts/mimic-tests/queries/inputs-cv.sql')
) %>% 
  select(
    -c(subject_id, hadm_id, itemid),
    charttime = charttime_equiv,
    value = amount,
    origin = orig_table
  ) %>% 
  as_tibble()

inputs_mv <- dbGetQuery(
  mimic,
  statement = read_file('scripts/mimic-tests/queries/inputs-mv.sql')
) %>% 
  select(
    -c(subject_id, hadm_id, itemid),
    charttime = charttime_equiv,
    value = amount,
    origin = orig_table
  ) %>% 
  as_tibble()

# outputs - this is conveniently all 'bolus' doses, I guess you don't really
# have fixed rate removal of fluids from a patient.
outputs <- dbGetQuery(
  mimic,
  statement = read_file('scripts/mimic-tests/queries/outputs.sql')
) %>%
  mutate(
    origin = 'outputs',
    value = -1 * value # make the outputs negative, so that we can average them all.
  ) %>%
  filter(value != 0) %>% 
  as_tibble()

# need to be since admission time, then we can do a window from there.
relevant_fluid_events <- bind_rows(
  inputs_cv,
  inputs_mv,
  outputs
) %>% 
  left_join(
    icustays %>% 
      select(-subject_id, hadm_id), 
    by = 'icustay_id'
  ) %>% 
  mutate(
    time_since_icu_adm = time_length(
      charttime - intime,
      unit = 'day'
    )
  ) %>% 
  as_tibble()
  
cumulative_fluids <- relevant_fluid_events %>%
  group_by(icustay_id) %>%
  arrange(icustay_id, time_since_icu_adm) %>%
  mutate(
    cumulative_fluids = cumsum(value)
  ) %>%
  select(-c(charttime, intime))


