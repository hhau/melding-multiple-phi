library(RPostgres)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)

source("scripts/common/logger-setup.R")
source("scripts/common/setup-argparse.R")

parser$add_argument("--combined-pf-and-summarised-fluid-data")
parser$add_argument("--demographics-query")
parser$add_argument("--median-labs-query")
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

combined_data <- readRDS(args$combined_pf_and_summarised_fluid_data)
icustay_ids <- combined_data$icustay_id %>%
  unique()

demographics <- dbGetQuery(
  mimic,
  read_file(args$demographics_query)
) %>%
  filter(icustay_id %in% icustay_ids)

first_labs <- dbGetQuery(
  mimic,
  read_file(args$median_labs_query)
) %>%
  filter(icustay_id %in% icustay_ids)

no_missing_values_subset_vec <- first_labs[, -(1 : 3)] %>%
  apply(2, function(x) sum(is.na(x))) == 0

first_labs_no_missing <- first_labs[, c(rep(TRUE, 3), no_missing_values_subset_vec)] %>%
  arrange(icustay_id)

baseline_covariate <- first_labs_no_missing %>%
  left_join(demographics, by = 'icustay_id')

saveRDS(
  file = args$output,
  object = baseline_covariate
)
