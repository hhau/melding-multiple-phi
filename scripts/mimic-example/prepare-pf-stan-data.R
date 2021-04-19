library(splines2)
library(dplyr)
library(abind)

source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")

parser$add_argument("--combined-pf-and-summarised-fluid-data")
parser$add_argument("--mimic-globals")
args <- parser$parse_args()

source(args$mimic_globals)

pf_data <- readRDS(args$combined_pf_and_summarised_fluid_data) %>%
  filter(value_type == 'pf')

split_pf_data <- pf_data %>%
  group_by(icustay_id) %>%
  arrange(icustay_id, time_since_icu_adm) %>%
  ungroup() %>%
  group_split(icustay_id)

res <- lapply(split_pf_data, function(a_patients_data) {
  boundary_knots <- c(
    min(a_patients_data$time_since_icu_adm),
    max(a_patients_data$time_since_icu_adm)
  )

  internal_knots <- seq(
    from = boundary_knots[1],
    to = boundary_knots[2],
    length.out = N_INTERNAL_KNOTS + 2
  )[-c(1, N_INTERNAL_KNOTS + 2)]

  x_obs_mat <- bSpline(
    x = a_patients_data$time_since_icu_adm,
    degree = 3,
    Boundary.knots = boundary_knots,
    knots = internal_knots,
    intercept = FALSE
  )

  x_plot_seq <- seq(
    from = boundary_knots[1],
    to = boundary_knots[2],
    length.out = N_PLOT_POINTS
  )

  x_plot_mat <- bSpline(
    x = x_plot_seq,
    degree = 3,
    Boundary.knots = boundary_knots,
    knots = internal_knots,
    intercept = FALSE
  )

  y_obs_vec <- a_patients_data$value
  y_obs_mean <- mean(y_obs_vec)
  y_obs_sd <- sd(y_obs_vec)
  y_obs_vec_ctr <- (y_obs_vec - y_obs_mean) / y_obs_sd
  threshold_ctr <- (GLOBAL_ARDS_THRESHOLD - y_obs_mean) / y_obs_sd

  res <- list(
    icustay_id = unique(a_patients_data$icustay_id),
    boundary_knots = boundary_knots,
    internal_knots = internal_knots,
    x_obs_mat = x_obs_mat,
    y_obs_vec = y_obs_vec,
    y_obs_mean = y_obs_mean,
    y_obs_sd = y_obs_sd,
    y_obs_vec_ctr = y_obs_vec_ctr,
    threshold_ctr = threshold_ctr,
    x_plot_seq = x_plot_seq,
    x_plot_mat = x_plot_mat
  )

  return(res)
})

bind_named_sublists <- function(outer_list, name, along_dim) {
  lapply(outer_list, function(sub_list) sub_list[[name]]) %>%
    abind(along = along_dim)
}

n_icu_stays <- length(res)

# each element gives the starting location of each patients data,
# except for the last element, which is the total number of rows + 1
subset_vector <- lapply(res, function(x) nrow(x$x_obs_mat)) %>%
  cumsum() %>%
  `+`(., 1) %>%
  c(1, .)

y_vec <- bind_named_sublists(res, "y_obs_vec", 1)
y_vec_ctr <- bind_named_sublists(res, "y_obs_vec_ctr", 1)
ctr_means <- bind_named_sublists(res, 'y_obs_mean', 1)
ctr_sds <- bind_named_sublists(res, "y_obs_sd", 1)
obs_matrix <- bind_named_sublists(res, "x_obs_mat", 1)
plot_matrix <- bind_named_sublists(res, "x_plot_mat", 3) %>%
  aperm(perm = c(3, 1, 2))

ids <- bind_named_sublists(res, "icustay_id", 1) %>%
  tibble(
    icustay_id = .,
    i = 1 : n_icu_stays
  )

res <- lapply(res, function(x) {
  x$i <- ids %>%
    filter(icustay_id == x$icustay_id) %>%
    pull(i)

  return(x)
})

n_basis_coef <- ncol(obs_matrix)
n_total_obs <- nrow(obs_matrix)

stan_data <- list(
  n_icu_stays = n_icu_stays,
  n_total_obs = n_total_obs,
  n_basis_coef = n_basis_coef,
  n_plot_points = N_PLOT_POINTS,
  subset_vector = subset_vector,
  obs_matrix = obs_matrix,
  y_vec_ctr = y_vec_ctr,
  ctr_means = ctr_means,
  ctr_sds = ctr_sds,
  plot_matrix = plot_matrix,
  spline_coef_prior_sd = 0.5,
  noise_df = 5
)

saveRDS(
  file = args$output,
  object = stan_data
)

list_format_name <- stringr::str_replace(
  args$output,
  'stan',
  'list'
)

saveRDS(
  file = list_format_name,
  object = res
)
