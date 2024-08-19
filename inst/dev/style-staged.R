staged_files <- system("git diff --cached --name-only --diff-filter=ACM", intern = TRUE)
if (length(staged_files) > 0) {
  staged_files <- staged_files[stringr::str_ends(string = staged_files, pattern = "\\.R")]
  if (length(staged_files) == 0) {
    cli::cli_alert_info("Found no staged R scripts to style.")
    quit(status = 0, save = "no")
  }
  res <- styler::style_file(path = staged_files)
  if (any(res$changed)) {
    cli::cat_rule()
    cat("\n")
    cli::cli_alert_warning(cli::col_br_blue("Stage the changed files by running:"))
    cat("git add", paste0(res$file[res$changed], collapse = " "), "\n")
  } else {
    cli::cli_alert_success("No changes made.")
  }
} else {
  cli::cli_alert_info("Found no staged R scripts to style")
  quit(status = 0, save = "no")
}
