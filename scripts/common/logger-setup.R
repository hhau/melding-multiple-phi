library(futile.logger)
library(crayon)

base_filename <- unlist(strsplit(x = basename(Sys.glob("./*.rmd")), ".rmd"))
log_filename <- paste0("logs/" , Sys.Date(), "_", base_filename, ".log") 
flog.appender(appender.file(log_filename), name = base_filename)

cat(
  blue(sprintf('futile.logger info')),
  "\n",
  "  ",
  blue$underline("log_filename:"), 
  sprintf("%s", log_filename),
  "\n",
  "  ",
  blue$underline("base_filename:"), 
  sprintf("%s", base_filename),
  "\n"
)

logfile_line_delineator <- function() {
  flog.info(
    paste0(paste0(rep("=", times = 40), collapse = ""), "\n"), 
    name = base_filename
  )
}