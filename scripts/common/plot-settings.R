library(ggplot2)
library(bayesplot)

# Theme settings
theme_set(theme_classic())
theme_replace(
  panel.grid.major = element_line(),
  panel.grid.minor = element_line(linetype = "dashed", size = rel(2/3)),
  legend.text = element_text(size = rel(1.1)),
  legend.title = element_text(size = rel(1.1)),
  strip.background = element_rect(fill = "#dee1e0"),
  strip.text = element_text(size = rel(1.0))
)

bayesplot_theme_set(theme_classic())
bayesplot_theme_replace(
  panel.grid.major = element_line(),
  panel.grid.minor = element_line(linetype = "dashed", size = rel(2/3)),
  legend.text = element_text(size = rel(1.1)),
  legend.title = element_text(size = rel(1.1)),
  strip.background = element_rect(fill = "#dee1e0"),
)

# Colours - should all be darkest[1] to lightest[n] 
blues <- c(
  "#00214F",
  "#2C7FB8",
  "#A6E9FF"
) 

greens <- c(
  "#364723",
  "#557438",
  "#749C4D",
  "#98CD65"
)

## Burgundy(ish) red highlight
highlight_col <- wesanderson::wes_palette("FantasticFox1")[5]

whw_pal_16 <- c(
  "#eab27e",
  "#f0ac84",
  "#f4a78c",
  "#f6a395",
  "#f59f9f",
  "#f19da9",
  "#ec9cb4",
  "#e49cbd",
  "#da9cc6",
  "#ce9ecd",
  "#c1a0d3",
  "#b2a2d7",
  "#a3a4d8",
  "#93a5d8",
  "#85a7d6",
  "#77a8d2"
)

# ggplot saving settings - @mbertolacci
display_settings <- list(
  full_page_plot_width = 15,
  full_page_plot_height = 21,
  half_page_plot_width = 7,
  half_page_plot_height = 10,
  png_plot_dpi = 300,
  highlight_colour = highlight_col
)

ggsave_base <- function(filename, plot, bg = 'transparent', ...) {
  ggsave(
    filename,
    plot,
    units = 'cm',
    dpi = display_settings$png_plot_dpi,
    bg = bg,
    ...
  )
}

ggsave_fullwidth <- function(filename, plot, ...) {
  ggsave_base(
    filename,
    plot,
    width = display_settings$full_page_plot_width,
    ...
  )
}

ggsave_fullpage <- function(filename, plot, adjust_height = 0, ...) {
  ggsave_base(
    filename,
    plot,
    width = display_settings$full_page_plot_width,
    height = display_settings$full_page_plot_height + adjust_height,
    ...
  )
}

ggsave_halfheight <- function(filename, plot, ...) {
  ggsave_base(
    filename,
    plot,
    width = display_settings$full_page_plot_width,
    height = display_settings$half_page_plot_height,
    ...
  )
}

parsed_map <- function(map) {
  function(x) parse(text = map[x])
}
