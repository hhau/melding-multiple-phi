library(RPostgres)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)

source("scripts/common/logger-setup.R")
source("scripts/common/setup-argparse.R")

parser$add_argument("--blood-gasses-query")
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

flog.info(
  'mimic-example: querying blood gasses',
  name = base_filename
)

res <- dbGetQuery(
  mimic,
  statement = read_file(args$blood_gasses_query)
)

tbl1 <- res %>%
  select(
    subject_id,
    hadm_id,
    icustay_id,
    charttime,
    specimen_pred,
    pf = pao2fio2
  ) %>%
  as_tibble() %>%
  filter(!is.na(pf), specimen_pred == 'ART') %>%
  drop_na()

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
)

# require the first n observations are above the threshold, avoids
# the case where the first or second observation is an outlier / clearly wrong
first_n_greater_than_k <- function(x, n, k) {
  all(x[1 : n] > k)
}

time_between_measurements_less_than <- function(x, gap_days = 2) {
  max(diff(x)) < gap_days
}

tbl2 <- tbl1 %>%
  left_join(icustays) %>%
  group_by(icustay_id) %>%
  mutate(
    time_since_icu_adm = time_length(
      charttime - intime,
      unit = 'day'
    )
  )

# initial filter speed up later computations? Possibly no longer needed
# now that we filter to at least 15 measurements.
min_obs <- 5
cohort_with_min_obs <- tbl2 %>%
  group_by(icustay_id) %>%
  count() %>%
  filter(n >= min_obs)

tbl3 <- tbl2 %>%
  filter(icustay_id %in% cohort_with_min_obs$icustay_id)

tbl4 <- tbl3 %>%
  group_by(icustay_id) %>%
  filter(
    pf < 600, ## infeasible pf ratios
    time_between_measurements_less_than(time_since_icu_adm, gap_days = 2),
    first_n_greater_than_k(pf, n = 6, k = 350)
  ) %>%
  count(icustay_id, sort = TRUE) %>%
  ungroup() %>%
  filter(between(n, 12, 500))

pf_tbl <- tbl2 %>%
  filter(icustay_id %in% tbl4$icustay_id)

flog.info(
  'mimic-example: writing pf data',
  name = base_filename
)

saveRDS(
  file = args$output,
  object = pf_tbl
)
