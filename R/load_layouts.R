#' Load layouts from an PXL file
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Layouts can be pre-computed with the Pixelator data processing pipeline and
#' are stored in a hive-styled partitioning in the PXL file. This function reads
#' the layouts from the PXL file and returns them as a list. Use
#' \code{\link{inspect_pxl_file}} to check the contents of a PXL file.
#'
#' @param filename Path to a PXL file
#' @param cells A character vector with component IDs. If NULL, all components are loaded.
#' @param graph_projection The graph projection to load. Default is 'bipartite'. If multiple
#'  projections are present in the file, only the selected one is loaded.
#' @param verbose Print messages
#'
#' @return A list of lists with the layouts. At the top level, the list is split by
#' layout. At the second level, the list is split by component. The components are sorted in
#' the order they appear in the PXL file.
#'
#' @export
#'
ReadMPX_layouts <- function(
  filename,
  cells = NULL,
  graph_projection = c("bipartite", "Anode", "linegraph"),
  verbose = TRUE
) {
  graph_projection <- match.arg(graph_projection,
    choices = c("bipartite", "Anode", "linegraph")
  )

  # Check file
  pxl_file_info <- inspect_pxl_file(filename)

  # Check cells
  stopifnot(
    "'cells' must be a non-empty character vector with component IDs" =
      is.null(cells) || (is.character(cells) && length(cells) > 0)
  )

  # Check if the file contains layouts
  if (!"layouts.parquet" %in% pxl_file_info$file_type) {
    abort(glue("File '{col_br_blue(filename)}' does not contain any layouts."))
  }

  # Create temporary directory
  temp_layout_dir <- .create_unique_temp_dir("temp_layouts")

  # Unzip layout parquet files
  layout_files <- pxl_file_info$file[[which(pxl_file_info$file_type == "layouts.parquet")]]
  cells <- cells %||% layout_files$component
  if (!all(cells %in% layout_files$component)) {
    abort(glue("File '{col_br_blue(filename)}' does not contain layouts for all components."))
  }
  layout_files <- layout_files %>%
    filter(component %in% cells)
  zip::unzip(filename, files = layout_files$file, exdir = temp_layout_dir)

  # Load hive-styled parquet files
  coords <- arrow::open_dataset(temp_layout_dir)
  coords <- coords %>%
    select(any_of(c("name", "x", "y", "z", "sample")), component, graph_projection, layout) %>%
    collect()

  # Handle graph_projections
  coords_layout_grouped <- coords %>%
    group_by(graph_projection)
  coords_layout_grouped_keys <- coords_layout_grouped %>%
    group_keys()

  # Keep selection graph_projection
  if (!graph_projection %in% coords_layout_grouped_keys$graph_projection) {
    abort(glue("Graph projection '{col_br_blue(graph_projection)}' does not exist in the file."))
  }
  if (!all(coords_layout_grouped_keys$graph_projection %in% graph_projection)) {
    coords <- coords %>%
      filter(graph_projection %in% graph_projection)
  }
  coords <- coords %>%
    select(-graph_projection)

  # Split coordinates by layout
  coords_layout_grouped <- coords %>%
    group_by(layout)
  coords_layout_split <- coords_layout_grouped %>%
    group_split() %>%
    set_names(nm = group_keys(coords_layout_grouped)$layout)

  # Split coordinates by component
  coords_layout_split <- lapply(coords_layout_split, function(coords_layout) {
    coords_component_grouped <- coords_layout %>%
      select(-layout) %>%
      group_by(component)
    coords_component_split <- coords_component_grouped %>%
      group_split() %>%
      lapply(function(x) {
        x %>%
          select(-component) %>%
          select(
            where(
              ~ sum(!is.na(.x)) > 0
            )
          )
      }) %>%
      set_names(nm = group_keys(coords_component_grouped)$component)
    return(coords_component_split[cells])
  })

  # Remove temporary directory
  err <- try(fs::dir_delete(temp_layout_dir))
  if (inherits(err, "try-error")) {
    warn(glue("Failed to remove temporary directory '{col_br_blue(temp_layout_dir)}'."))
  }

  return(coords_layout_split)
}
