library(tidyverse)
library(Rfast)

n_data <- 5e5

# the code assumes the times are numeric. If they are dates,
# the coerce them to numeric by lubridate::time_length,
# with the appropriate start time (usually admission), end time
# and unit (usually days, sometimes hours for shorter stays).
times <- c(runif(n = n_data, min = -1.5, max = 3)) %>%
  sort()

values <- arima.sim(
  model = list(ar = c(0.5, 0.2), ma = c(0.2)),
  n = n_data
)

df <- tibble(times, values)

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

plot_origin <- 0
plot_interval <- 6/24

plot_df <- df %>% mutate(
  grp = make_window_groups(times, plot_origin, plot_interval) %>%
    as.factor()
)

boundary_df <- df %>%
  pull(times) %>%
  make_window_groups(plot_origin, plot_interval, return_boundaries = TRUE) %>%
  tibble(value = .)

ggplot(plot_df, aes(x = times, y = values, col = grp)) +
  geom_point() +
  geom_vline(
    data = boundary_df,
    mapping = aes(xintercept = value)
  ) +
  ylab(expression(Delta ~ 'fluid')) +
  xlab('time')

summarised_df <- plot_df %>%
  group_by(grp) %>%
  summarise(
    values = sum(values),
    times = mean(times) # where should the time go? Not hugely important.
  )

ggplot(summarised_df, aes(x = times, y = values, col = grp)) +
  geom_point() +
  geom_vline(data = boundary_df, mapping = aes(xintercept = value)) +
  ylab(expression(Delta ~ 'fluid')) +
  xlab('time')

