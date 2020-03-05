library(dplyr)
library(tibble)
library(purrr)

base_dir <- "other-code/monte-carlo-rcpp/projects/recapture/examples/owls/data"

# Fecundity
fecundity_tbl <- read.table(
  file = file.path(base_dir, "fecundity.dat")
) %>% 
  as_tibble() %>% 
  transmute(
    N_breeding_females = V1,
    N_offspring = V2
  )

saveRDS(
  object = fecundity_tbl,
  file = "rds/owls-example/fecundity-data.rds"
)

# Count data
count_tbl <- read.table(
  file = file.path(base_dir, "count.dat")
) %>% 
  as_tibble() %>% 
  transmute(y_count_female = V1)

saveRDS(
  object = count_tbl,
  file = "rds/owls-example/count-data.rds"
)

# Capture-Recapture
base_name <- "capRecap"
gender <- c("Female", "Male")
age <- c("First", "Adult")

invisible(lapply(gender, function(.x) {
  lapply(age, function(.y) {
    data_file_name <- paste0(base_name, .x, .y, ".dat")
    file_path <- file.path(base_dir, data_file_name)
    a_matrix <- read.table(file_path) %>% 
      as.matrix()  
      
    colnames(a_matrix) <- NULL
    
    saveRDS(
      object = a_matrix,
      file = file.path(
        "rds/owls-example",
        paste0(
          paste(
            "capture-recapture",
            tolower(.x),
            tolower(.y),
            sep = "-"
          ),
        "-data.rds")
      )
    )
  })
}))
