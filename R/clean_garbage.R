# Declarations used in package check
globalVariables(
  names = c('edgelist_dir', 'variable', 'pattern_match', 'type', 'size'),
  package = 'pixelatorR',
  add = TRUE
)

#' Clean up Unlinked Edgelist Directories
#'
#' Each time a \code{Seurat} with a \code{CellGraphAssay} or a \code{CellGraphAssay} is created,
#' renamed, filtered (by cells) or merged, a copy of the associated edgelist is made
#' in the \code{getOption("pixelatorR.arrow_outdir")} directory. For large data sets,
#' these edgelist files can use up several gigabytes of disk space. This function
#' scans through the \code{getOption("pixelatorR.arrow_outdir")} directory and removes
#' any edgelist files that are not linked to \code{Seurat} objects or \code{CellGraphAssay}s
#' in the global environment. Note that only global variables which inherits the \code{Seurat}
#' or \code{CellGraphAssay} class are considered, as well as lists containing these object types.
#'
#' @export
#'
clean_edgelists_directories <- function () {

  # Get variables in global environment
  global_variables <- global_env() %>% as.list()

  # Check if there are any Seurat objects or CellGraphAssays
  global_variable_types <- sapply(names(global_variables), function(nm) {

    x <- global_variables[[nm]]

    if (inherits(x, what = c("Seurat", "CellGraphAssay"))) {
      # Ignore Seurat objects without CellGraphAssays
      arrow_dir <- try({
        ArrowDir(x)
      }, silent = TRUE)
      if (inherits(arrow_dir, "try-error")) {
        return(NULL)
      } else {
        return(arrow_dir)
      }
    } else if (inherits(x, what = "list")) {

      # Check if any list contains Seurat objects or CellGraphAssays
      list_arrow_dirs <- lapply(x, function(y) {
        if (inherits(y, what = c("Seurat", "CellGraphAssay"))) {
          # Ignore Seurat objects without CellGraphAssays
          arrow_dir <- try({
            ArrowDir(y)
          }, silent = TRUE)
          if (inherits(arrow_dir, "try-error")) {
            return(NULL)
          } else {
            return(arrow_dir)
          }
        } else {
          return(NULL)
        }
      })
      return(list_arrow_dirs %>% set_names(seq_along(list_arrow_dirs)))
    } else {
      return(NULL)
    }
  })

  # Unlist and remove NULLs
  global_variable_types <- unlist(global_variable_types)
  global_variable_types <- Filter(Negate(is.null), global_variable_types)

  # Convert to tibble
  if (length(global_variable_types) > 0) {
    global_variable_types <- tibble(variable = names(global_variable_types), edgelist_dir = global_variable_types)
  }

  # Get command log
  command_log <- getOption("pixelatorR.edgelist_copies")
  if (length(command_log) > 0) {
    if (length(global_variable_types) > 0) {
      command_log <- command_log %>%
        left_join(global_variable_types, by = "edgelist_dir")
    } else {
      command_log <- NULL
    }
  }

  # Get edgelist directories
  edgelist_dirs <- list.files(getOption("pixelatorR.arrow_outdir"), full.names = TRUE) %>% normalizePath()
  if (length(edgelist_dirs) > 0) {
    edgelist_dirs <- tibble(edgelist_dir = edgelist_dirs, exists = TRUE)
    if (length(global_variable_types) > 0) {
      command_log <- command_log %>%
        full_join(edgelist_dirs, by = "edgelist_dir")
    } else {
      command_log <- edgelist_dirs
    }
  } else {
    cli_alert_info("Found no edgelist directories linked to missing variables")
    return(invisible(NULL))
  }

  # Only keep existing directories
  command_log <- command_log %>%
    filter(dir.exists(edgelist_dir))

  if (!"variable" %in% names(command_log)) {
    command_log <- command_log %>% mutate(variable = NA_character_)
  }

  # Remove edgelists linked to missing variables
  command_log_na <- command_log %>%
    filter(is.na(variable))

  # Run check on directories to make sure that they match the expected pattern
  command_log_na <- command_log_na %>%
    mutate(
      pattern_match = stringr::str_detect(basename(edgelist_dir),
                                          pattern = "[A-Za-z0-9]{5}-\\d{4}-\\d{2}-\\d{2}-\\d{6}")
    ) %>%
    filter(pattern_match)

  # Exit if no directories are found
  if (nrow(command_log_na) == 0) {
    cli_alert_info("Found no edgelist directories linked to missing variables")
    return(invisible(NULL))
  }

  if (getOption("pixelatorR.interactive", default = TRUE)) {
    cli_alert_warning(glue("Are you sure you want to remove the following directories?\n\n",
                           paste(command_log_na$edgelist_dir, collapse = "\n")))
    if (menu(c("Yes", "No")) != 1) {
      return(invisible(NULL))
    }
  }

  file_sizes <- c()

  # Remove directories
  cli_alert_info(glue("Removing {nrow(command_log_na)} edgelist directories linked to missing variables"))
  for (d in unique(command_log_na$edgelist_dir)) {
    file_sizes <- c(file_sizes, fs::dir_info(d, recurse = TRUE) %>%
                      filter(type == "file") %>%
                      pull(size))
    fs::dir_delete(path = d)
  }

  cli_alert_info("Freed up {sum(file_sizes) %>% fs::fs_bytes()} disk space")

  # Update global option pixelatorR.edgelist_copies
  command_log_clean <- command_log %>%
    filter(!is.na(variable)) %>%
    select(-variable, -exists)
  options(pixelatorR.edgelist_copies = command_log_clean)

}
