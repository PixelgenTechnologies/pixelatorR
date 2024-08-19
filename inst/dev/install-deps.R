options(repos = c(CRAN = "https://cloud.r-project.org"))
if (getRversion() < "4.4.1") {
  cli_alert_error("Please upgrade to R 4.4.1 or later.")
  quit(status = 1)
}

if (!requireNamespace("pak", quietly = TRUE)) {
  cat("Installing pak.\n")
  cat("This can take a few minutes...\n")
  install.packages("pak", quiet = TRUE)
}
if (!requireNamespace("desc", quietly = TRUE)) {
  pak::pkg_install(c("desc", "cli"))
}
cli::cli_h3("Attempting to install hdf5r...")
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  pak::pkg_install("cran/hdf5r")
}
cli::cli_h2("Installing pixelatorR dependencies")
find_deps <- desc::description$new("DESCRIPTION")$get_deps()
pak::pkg_install(pkg = find_deps$package, ask = FALSE)
cli::cli_h2("Installing pixelatorR")
pak::pkg_install(".", ask = FALSE, dependencies = FALSE)
cli::cat_line()
cli::cli_rule()
cli::cli_alert_success("Successfully installed pixelatorR!")
