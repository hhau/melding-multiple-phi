library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(pbapply)

source("scripts/common/plot-settings.R")

phi_12_samples <- readRDS("rds/surv-example/stage-one-phi-12-samples.rds")
n_iter <- dim(phi_12_samples)[1]
event_time_names <- names(phi_12_samples[1, 1,]) %>%
  grep("time", x = ., value = T)

event_indicator_names <- names(phi_12_samples[1, 1, ]) %>%
  grep("indicator", x = ., value = T)

# it would be nice to speed this up some, or at least cache the results
res <- pblapply(1 : n_iter, cl = 6, function(iter_id) {
  df <- data.frame(
    event_time = phi_12_samples[iter_id, 1, event_time_names],
    event_indicator = phi_12_samples[iter_id, 1, event_indicator_names],
    iter_id = iter_id
  )
  surv_obj <- survfit(Surv(event_time, event_indicator) ~ 1, data = df) %>% 
    surv_summary() %>%
    as_tibble() %>%
    select(-c(std.err, upper, lower)) %>%
    mutate(iter_id = iter_id)
  
  base_df <- surv_obj[1, ] %>%
    mutate(time = 0, surv = 1, n.event = 0)
  
  res <- rbind(base_df, surv_obj)
}) 

plot_df <- bind_rows(res)
censor_df <- plot_df %>% 
  filter(n.censor > 0) %>%
  mutate(plot_alpha = (n.censor - min(n.censor)) / (diff(range(n.censor))))

random_draws <- sample(1 : n_iter, size = 2000)

ggplot(
  data = plot_df %>% filter(iter_id %in% random_draws), 
  aes(x = time, y = surv, group = iter_id)
) +
  geom_step(alpha = 0.5) + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(
    data = censor_df %>% filter(iter_id %in% random_draws),
    inherit.aes = FALSE,
    mapping = aes(x = time, y = surv),
    alpha = censor_df %>% 
      filter(iter_id %in% random_draws) %>%
      pull(plot_alpha),
    shape = 3,
    col = "red"
  )

# add hazard based surv prob -- note that this isn't correct
# yet because it doesn't account for variation from the linear predictor / 
# scale parameter. If we want to include this then I need to think about it
# A visually rough estimate could be had by considering the surv_prob_term_one
# and surv_prob_term_two, associated with each hazard draw
psi_2_samples <- readRDS(file = "rds/surv-example/stage-two-psi-2-samples.rds")
haz_samples <- psi_2_samples[, , "hazard_gamma"] %>% as.numeric()
plot_t <- seq(from = 0, to = 1, length.out = 25)
haz_df <- expand.grid(
  plot_t = plot_t,
  hazard_gamma = haz_samples
)

haz_df$surv_prob <- pweibull(haz_df$plot_t, shape = haz_df$hazard_gamma, lower.tail = FALSE)

haz_plot_tbl <- haz_df %>%
  as_tibble() %>%
  select(-hazard_gamma) %>%
  group_by(plot_t) %>%
  summarise(
    mean = mean(surv_prob),
    median = median(surv_prob),
    lower = quantile(surv_prob, 0.1),
    upper = quantile(surv_prob, 0.9)
  ) %>%
  mutate(method = "chained")

ggplot(haz_plot_tbl, aes(x = plot_t)) +
  geom_line(aes(y = mean), col = "blue") +
  geom_line(aes(y = median), col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

all_phi_fixed_samples <- readRDS("rds/surv-example/point-est-psi-2-samples.rds")
fixed_haz_samples <- all_phi_fixed_samples[, , "hazard_gamma"] %>% as.numeric()
fixed_haz_df <- expand.grid(
  plot_t = plot_t,
  hazard_gamma = fixed_haz_samples
)

fixed_haz_df$surv_prob <- pweibull(
  fixed_haz_df$plot_t, 
  shape = fixed_haz_df$hazard_gamma, 
  lower.tail = FALSE
)

fixed_haz_plot_tbl <- fixed_haz_df %>%
  as_tibble() %>%
  select(-hazard_gamma) %>%
  group_by(plot_t) %>%
  summarise(
    mean = mean(surv_prob),
    median = median(surv_prob),
    lower = quantile(surv_prob, 0.1),
    upper = quantile(surv_prob, 0.9)
  ) %>%
  mutate(method = "all fixed")

compare_haz_tbl <- rbind(haz_plot_tbl, fixed_haz_plot_tbl) %>%
  mutate(method = as.factor(method))

ggplot(compare_haz_tbl, aes(x = plot_t, col = method, group = method)) +
  geom_line(aes(y = mean)) +
  # geom_line(aes(y = median), col = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1))

# gross big hacky plot:
ggplot(
  data = plot_df %>% filter(iter_id %in% random_draws), 
  aes(x = time, y = surv, group = iter_id)
) +
  geom_step(alpha = 0.025) + 
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(
    data = censor_df %>% filter(iter_id %in% random_draws),
    inherit.aes = FALSE,
    mapping = aes(x = time, y = surv),
    alpha = censor_df %>% 
      filter(iter_id %in% random_draws) %>%
      pull(plot_alpha),
    shape = 3,
    col = "red"
  ) +
  geom_line(
    inherit.aes = FALSE,
    data = compare_haz_tbl,
    aes(x = plot_t, y = mean, col = method, group = method), 
  ) +
  # geom_line(
  #   inherit.aes = FALSE,
  #   data = compare_haz_tbl,
  #   aes(x = plot_t, y = median, col = method, group = method), 
  #   col = "blue",
  #   linetype = "dashed"
  # ) +
  geom_ribbon(
    inherit.aes = FALSE,
    data = compare_haz_tbl,
    aes(x = plot_t, ymin = lower, ymax = upper, col = method, group = method), 
    alpha = 0.2,
  )

