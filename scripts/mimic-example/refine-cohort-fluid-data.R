library(dplyr)
library(Rfast)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--combined-pf-and-raw-fluid-data")
parser$add_argument("--mimic-globals")
args <- parser$parse_args()

source(args$mimic_globals)

pf_and_raw_fluid_data <- readRDS(args$combined_pf_and_raw_fluid_data)
raw_fluid_data <- pf_and_raw_fluid_data %>%
  filter(value_type == 'fluids')

pf_data <- pf_and_raw_fluid_data %>%
  filter(value_type == 'pf')

make_window_groups <- function(
  t,
  origin = 0,
  interval = 1,
  return_boundaries = FALSE
) {
  n_t <- length(t)
  min_t <- min(t)
  max_t <- max(t)

  lower_boundary <- floor(min_t - interval)
  upper_boundary <- ceiling(max_t + interval)

  sorted_group_boundaries <- c(
    seq(from = origin, to = upper_boundary, by = interval),
    seq(from = origin, to = lower_boundary, by = -interval)
  ) %>%
    sort_unique() %>%
    subset(., between(., min_t, max_t))

  if (return_boundaries) {
    return(sorted_group_boundaries)
  }

  n_boundaries <- length(sorted_group_boundaries)
  grouping_list <- list()

  for (ii in 1 : n_boundaries) {
    group_vec <- array(data = n_boundaries + 1, dim = n_t)
    indicies <- which(t < sorted_group_boundaries[ii])
    group_vec[indicies] <- ii
    grouping_list[[ii]] <- group_vec
  }

  res <- grouping_list %>%
    simplify2array() %>%
    apply(1, min)

  return(res)

}

summarised_fluid_data <- raw_fluid_data %>%
  group_by(icustay_id) %>%
  mutate(
    grp = make_window_groups(
      t = time_since_icu_adm,
      origin = FLUID_TIME_ORIGIN,
      interval = FLUID_WINDOW_WIDTH_DAYS
    ) %>%
      as.factor()
  ) %>%
  group_by(icustay_id, grp) %>%
  summarise(
    value = sum(value),
    time_since_icu_adm = mean(time_since_icu_adm)
  ) %>%
  select(-grp) %>%
  mutate(value_type = 'fluids')

final_data <- bind_rows(pf_data, summarised_fluid_data) %>%
  select(-c(amountuom, label))

saveRDS(
  file = args$output,
  object = final_data
)
