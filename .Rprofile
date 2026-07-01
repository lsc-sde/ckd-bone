source("renv/activate.R")
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  stringsAsFactors = FALSE
)

if (requireNamespace("renv", quietly = TRUE)) {
  tryCatch(
    renv::activate(),
    error = function(e) message("renv activate skipped: ", e$message)
  )
}
