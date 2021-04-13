library(dplyr)
library(ggplot2)
library(ggh4x)
library(mgcv)

source("scripts/common/plot-settings.R")
source("scripts/common/setup-argparse.R")
source("scripts/common/logger-setup.R")
source("scripts/mimic-example/GLOBALS.R")

parser$add_argument("--combined-pf-and-summarised-fluid-data")
args <- parser$parse_args()

plot_data <- readRDS(args$combined_pf_and_summarised_fluid_data)

psi_t <- function(x, nu = 3, deriv = 0) {
  if (deriv == 0) {
    res <- 1 / (x^2 + nu)
  } else if (deriv == 1) {
    res <- (nu - x^2) / ((x^2 + nu)^2)
  }

  return(res)
}

p1 <- ggplot(plot_data, aes(x = time_since_icu_adm, y = value)) +
  geom_point() +
  facet_nested_wrap(vars(icustay_id, value_type), scales = 'free_y') +
  stat_smooth(
    data = plot_data %>%
      filter(value_type == 'pf'),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    method.args = list(
      family = mgcv::scat(min.df = 5)
    )
  ) +
  geom_hline(
    data = tibble(
      value_type = 'pf',
      yintercept = 300
    ),
    mapping = aes(
      yintercept = yintercept
    )
  ) +
  stat_smooth(
    data = plot_data %>%
      filter(value_type == 'fluids'),
    method = MASS::rlm,
    method.args = list(psi = psi_t)
  ) +
  xlab('Days since start of this ICU admission') +
  ggtitle(
    label = sprintf(
      'PF ratio and net fluid change for a %d day window',
      FLUID_WINDOW_WIDTH_DAYS
    ),
    subtitle = sprintf(
      'Plot/data as of %s',
      Sys.time()
    )
  ) +
  theme(
    strip.text = element_text(margin = margin(t = 2, b = 2, l = 50, r = 50)),
    axis.title.y = element_blank()
  )

ggsave(
  filename = args$output,
  plot = p1,
  height = 30,
  width = 45,
  units = 'in'
)
