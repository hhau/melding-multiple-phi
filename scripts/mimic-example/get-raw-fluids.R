library(RPostgres)
library(dplyr)
library(readr)
library(lubridate)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--inputs-cv-query")
parser$add_argument("--inputs-mv-query")
parser$add_argument("--outputs-query")
args <- parser$parse_args()

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
      mimiciii.icustays;"
) %>%
  as_tibble()

flog.info(
  "mimic-example: querying carevue fluid data",
  name = base_filename
)

inputs_cv <- dbGetQuery(
  mimic,
  statement = read_file(args$inputs_cv_query)
) %>%
  select(
    -c(subject_id, hadm_id, itemid),
    charttime = charttime_equiv,
    value = amount,
    origin = orig_table
  ) %>%
  as_tibble()

flog.info(
  "mimic-example: querying metavision fluid data",
  name = base_filename
)

inputs_mv <- dbGetQuery(
  mimic,
  statement = read_file(args$inputs_mv_query)
) %>%
  select(
    -c(subject_id, hadm_id, itemid),
    charttime = charttime_equiv,
    value = amount,
    origin = orig_table
  ) %>%
  as_tibble()

flog.info(
  "mimic-example: querying outputs fluid data",
  name = base_filename
)

outputs <- dbGetQuery(
  mimic,
  statement = read_file(args$outputs_query)
) %>%
  mutate(
    value = -1 * value # make the outputs negative, so that we can average them all.
  ) %>%
  filter(value != 0) %>%
  rename(origin = orig_table) %>%
  as_tibble()

# Set time origin as admission time, then we can do a window from there.
# filter to complete cases as, as there are a handful of NA values
# due to bag/bottle changes with non-zero entries
# and missing admission identifiers.
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
  filter(complete.cases(.)) %>%
  as_tibble()

flog.info(
  "mimic-example: writing raw fluid data",
  name = base_filename
)

saveRDS(
  file = args$output,
  object = relevant_fluid_events
)
