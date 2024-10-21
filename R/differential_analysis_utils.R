#' Select columns from a \code{tbl_df} for grouping
#'
#' @returns A \code{tbl_df} object
#'
#' @noRd
.select_group_data <- function(object, contrast_column, group_vars = NULL) {
  group_data <- object %>%
    {
      if (!is.null(group_vars)) {
        stopifnot(
          "'group_vars' must be valid meta data column names" =
            inherits(group_vars, what = "character") &&
            (length(group_vars) >= 1) &&
            all(group_vars %in% colnames(.))
        )
        select(., all_of(c(contrast_column, group_vars)))
      } else {
        select(., all_of(contrast_column))
      }
    } %>%
    as_tibble(rownames = "component")
}

#' Prints the number of groups and comparisons for a differential analysis
#'
#' @noRd
.print_da_group_strategy <- function(test_groups_keys, targets, reference) {
  if (length(targets) > 1) {
    cli_alert_info(glue(
      "Running {nrow(test_groups_keys)} tests for each ",
      "of the following comparisons:\n"
    ))
  } else {
    cli_alert_info(glue(
      "Running {nrow(test_groups_keys)} tests for the following comparison:\n"
    ))
  }
  print(glue("  - {paste(targets, reference, sep = ' vs ')}"))
}
