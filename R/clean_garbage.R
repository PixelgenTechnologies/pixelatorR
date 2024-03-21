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

  global_variable_types <- .get_global_mpx_variables()

  # Convert to tibble
  if (length(global_variable_types) > 0) {
    global_variable_types <- tibble(variable = names(global_variable_types), edgelist_dir = global_variable_types)
  }

  # Get edgelist directories
  edgelist_dirs <- list.files(getOption("pixelatorR.arrow_outdir"), full.names = TRUE) %>% normalizePath()
  if (length(edgelist_dirs) > 0) {
    edgelist_dirs <- tibble(edgelist_dir = edgelist_dirs, exists = TRUE)
    if (length(global_variable_types) > 0) {
      global_variable_types <- global_variable_types %>%
        full_join(edgelist_dirs, by = "edgelist_dir")
    } else {
      global_variable_types <- edgelist_dirs
    }
  } else {
    cli_alert_info("Found no edgelist directories linked to missing variables")
    return(invisible(NULL))
  }

  # Only keep existing directories
  global_variable_types <- global_variable_types %>%
    filter(dir.exists(edgelist_dir))

  if (!"variable" %in% names(global_variable_types)) {
    global_variable_types <- global_variable_types %>% mutate(variable = NA_character_)
  }

  # Remove edgelists linked to missing variables
  global_variable_types_na <- global_variable_types %>%
    filter(is.na(variable))

  # Run check on directories to make sure that they match the expected pattern
  global_variable_types_na <- global_variable_types_na %>%
    mutate(
      pattern_match = stringr::str_detect(basename(edgelist_dir),
                                          pattern = "[A-Za-z0-9]{5}-\\d{4}-\\d{2}-\\d{2}-\\d{6}")
    ) %>%
    filter(pattern_match)

  # Exit if no directories are found
  if (nrow(global_variable_types_na) == 0) {
    cli_alert_info("Found no edgelist directories linked to missing variables")
    return(invisible(NULL))
  }

  if (getOption("pixelatorR.auto_cleanup", default = TRUE)) {
    cli_alert_warning(glue("Are you sure you want to remove the following directories?\n\n",
                           paste(global_variable_types_na$edgelist_dir, collapse = "\n")))
    if (menu(c("Yes", "No")) != 1) {
      return(invisible(NULL))
    }
  }

  file_sizes <- c()

  # Remove directories
  cli_alert_info(glue("Removing {nrow(global_variable_types_na)} edgelist ",
                      "directories linked to missing variables"))
  for (d in unique(global_variable_types_na$edgelist_dir)) {
    file_sizes <- c(file_sizes, fs::dir_info(d, recurse = TRUE) %>%
                      filter(type == "file") %>%
                      pull(size))
    if (fs::dir_exists(d)) {
      fs::dir_delete(path = d)
    }
  }

  cli_alert_info("Freed up {sum(file_sizes) %>% fs::fs_bytes()} disk space")

}


#' Check Disk Usage of Edgelist Directories
#'
#' This function checks the directories in \code{getOption("pixelatorR.arrow_outdir")}
#' and returns their total disk usage.
#'
#' @param list_all Show the disk usage for each sub directory and check what global variable,
#' if any, the sub directories are associated with.
#'
#' @return The total disk usage of all valid edgelist directories or a tibble
#' with the disk usage of each sub directory and the variable they are linked
#' to if applicable.
#'
#' @export
#'
edgelist_directories_du <- function (
  list_all = FALSE
) {

  # Fetch all edgelist directories
  edgelist_dirs <- list.files(getOption("pixelatorR.arrow_outdir"), full.names = TRUE) %>% normalizePath()
  edgelist_dirs <- edgelist_dirs[stringr::str_detect(basename(edgelist_dirs),
                      pattern = "[A-Za-z0-9]{5}-\\d{4}-\\d{2}-\\d{2}-\\d{6}")]

  # Create a tibble with the disk usage of each sub directory
  file_info <- lapply(edgelist_dirs, function(f) {
    fs::dir_info(f, recurse = TRUE) %>%
      filter(type == "file") %>%
      select(path, size) %>%
      mutate(path = path %>% dirname() %>% dirname())
  }) %>%
    do.call(bind_rows, .)
  if (list_all) {

    # Fetch global variables linked to edgelists
    global_mpx_variables <- .get_global_mpx_variables()
    global_mpx_variables <- tibble(variable = names(global_mpx_variables),
                                   path = global_mpx_variables)

    # If there are any global variables linked to edgelists, join them to the file_info tibble
    if (nrow(global_mpx_variables) > 0) {
      file_info <- left_join(file_info, global_mpx_variables, by = "path")
    }
    if (nrow(file_info) > 0) {
      return(file_info)
    } else {
      cli_alert_info("Found no edgelist directories.")
      return(invisible(NULL))
    }
  } else {
    # Default behavior is to return the total disk usage
    # Note that some directories might be duplicated if
    # they are linked to multiple variables
    return(fs::fs_bytes(sum(file_info %>% filter(!duplicated(path)) %>% pull(size))))
  }
}


#' Get Global Variables Linked to Edgelists
#'
#' @noRd
.get_global_mpx_variables <- function () {

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

  return(global_variable_types)
}
