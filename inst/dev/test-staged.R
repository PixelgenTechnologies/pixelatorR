staged_files <- system("git diff --cached --name-only --diff-filter=ACM", intern = TRUE)
if (length(staged_files) > 0) {
  staged_files <- staged_files[stringr::str_ends(string = staged_files, pattern = "\\.R")]
  staged_files <- staged_files[stringr::str_starts(string = basename(staged_files), pattern = "test-")]
  if (length(staged_files) == 0) {
    cli::cli_alert_info("Found no staged R test scripts")
    quit(status = 0, save = "no")
  }
  for (file in staged_files) {
    devtools::test_active_file(file = file)
  }
} else {
  cli::cli_alert_info("Found no staged files")
  quit(status = 0, save = "no")
}
