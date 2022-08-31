library(magrittr)

lower_case_string_to_integer <- function(x) {
  x %>% 
    strsplit(split = "") %>% 
    unlist() %>% 
    sapply(function(x) which(letters == x)) %>% 
    paste0(collapse = "") %>% 
    as.numeric() %>% 
    `%%`(.Machine$integer.max) %>% 
    as.integer()
}

data_seed <- "dataseed" %>% 
  lower_case_string_to_integer()

sim_seed <- "simseed" %>% 
  lower_case_string_to_integer()
  
# Global settings
n_patients <- 36L
n_long_beta <- 2L