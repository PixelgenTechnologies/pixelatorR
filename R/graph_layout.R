#' @include generics.R
NULL

#' Compute graph layouts
#'
#' Computes graph layouts for component graphs.
#'
#' @param layout_method The method for calculating the graph layout:
#' weighted Pivot MDS ("wpmds") or Pivot MDS (pmds)
#' @param dim An integer specifying the dimensions of the layout (2 or 3). Note that
#' for "pmds" and "wpmds", the x and y coordinates will be identical regardless if
#' dim is set to 2 or 3.
#' @param normalize_layout Logical specifying whether the coordinate system
#' should be centered at origo and the coordinates scaled such that their median
#' length (euclidean norm) is 1.
#' @param pivots Only used for "wpmds" and "pmds" layout algorithms.
#' See \code{?layout_with_pmds} for details
#' @param project_on_unit_sphere Should the resulting layout be projected onto
#'  a unit sphere?
#' @param seed Set seed for reproducibility
#' @param custom_layout_function A custom function for layout computation. The
#' function should take a \code{tbl_graph} object and a \code{dim} value as
#' input and return a matrix of dimensions NxD, where N is the number of nodes
#' and D is equal to \code{dim}. Note that this will override the \code{layout_method}.
#' @param custom_layout_function_args A list of arguments passed to \code{custom_layout_function}.
#' The \code{dim} is automatically passed to \code{custom_layout_function} and should not be
#' included in \code{custom_layout_function_args}.
#'
#' @rdname ComputeLayout
#' @method ComputeLayout tbl_graph
#'
#' @seealso [center_layout_coordinates()] for centering layout coordinates,
#' [normalize_layout_coordinates()] for normalizing layout coordinates and
#' [project_layout_coordinates_on_unit_sphere()] for projecting layout coordinates onto a unit sphere.
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_mpx_pxl_file()
#'
#' # Load example data
#' seur <- ReadMPX_Seurat(pxl_file)
#'
#' # Load 1 cellgraph
#' seur <- LoadCellGraphs(seur,
#'   cells = colnames(seur)[1],
#'   force = TRUE
#' )
#'
#' # Get CellGraph
#' cg <- CellGraphs(seur)[[colnames(seur)[1]]]
#'
#' # Get tbl_graph object
#' tbl_graph <- slot(cg, name = "cellgraph")
#'
#' # Compute layout for a tbl_graph
#' layout <- ComputeLayout(tbl_graph, layout_method = "wpmds")
#' layout %>% head()
#'
#' @export
#'
ComputeLayout.tbl_graph <- function(
  object,
  layout_method = c("pmds", "wpmds"),
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  ...
) {
  # Validate input parameters
  assert_single_value(dim, type = "numeric")
  assert_single_value(normalize_layout, type = "bool")
  assert_single_value(seed, type = "integer")

  if (project_on_unit_sphere && dim == 2) {
    cli::cli_abort(
      c(
        "i" = "Projecting onto a unit sphere is only possible for 3D layouts",
        "x" = "Change {.var dim} or set {.var project_on_unit_sphere=FALSE}"
      )
    )
  }
  if (project_on_unit_sphere && normalize_layout) {
    cli::cli_abort(
      c("x" = "Only one of {.var project_on_unit_sphere} or {.var normalize_layout} can be set to TRUE")
    )
  }

  # Set seed
  set.seed(seed)

  # validate and use custom layout function if available
  if (!is.null(custom_layout_function)) {
    .validate_custom_layout_function(custom_layout_function, custom_layout_function_args, dim)
    layout_function <- custom_layout_function
    layout_method <- "custom"

    # Try custom layout function on tbl_graph
    layout <- try({
      do.call(custom_layout_function, c(list(g = object), custom_layout_function_args))
    })
    if (inherits(layout, what = "try-error")) {
      cli::cli_abort(
        c("x" = "{.var custom_layout_function} failed to compute layout")
      )
    }

    .validate_custom_layout_function_results(layout, dim, n_nodes = length(object))
  } else {
    # Check and select a layout method
    layout_method <- match.arg(layout_method, choices = c("pmds", "wpmds"))

    layout_function <-
      switch(layout_method,
        "wpmds" = layout_with_weighted_pmds,
        "pmds" = {
          # Abort if pmds is selected and graphlayouts isn't installed
          expect_graphlayouts()
          graphlayouts::layout_with_pmds
        }
      )

    layout <-
      object %>%
      layout_function(dim = dim, pivots = pivots, ...) %>%
      as_tibble(.name_repair = function(x) c("x", "y", "z")[seq_len(dim)]) %>%
      bind_cols(object %N>% as_tibble()) %>%
      select(any_of(c("x", "y", "z")))
  }

  # Project to unit sphere
  if (project_on_unit_sphere) {
    layout <- project_layout_coordinates_on_unit_sphere(layout)
  }

  # Normalize layout
  if (normalize_layout) {
    layout <- normalize_layout_coordinates(layout)
  }

  return(layout)
}

#' @param layout_name The name of the computed layout. If this name is not given,
#' the \code{layout_method} will be used as the name. If \code{dim = 3}, a suffix
#' of "_3d" will be added to the layout name. This behavior is ignored if \code{layout_name}
#' is provided.
#'
#' @rdname ComputeLayout
#' @method ComputeLayout CellGraph
#'
#' @examples
#'
#' # Compute layout for a CellGraph
#' cg <- ComputeLayout(cg, layout_method = "wpmds")
#'
#' @export
#'
ComputeLayout.CellGraph <- function(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  ...
) {
  if (is.null(custom_layout_function) && is.null(layout_name)) {
    layout_name <- match.arg(layout_method, choices = c("pmds", "wpmds"))
    # Add suffix _3d if dim = 3
    if (dim == 3) {
      layout_name <- paste0(layout_name, "_3d")
    }
  }
  if (!is.null(custom_layout_function) && is.null(layout_name)) {
    layout_name <- "custom"
  }

  assert_single_value(layout_name, type = "string", allow_null = TRUE)

  layout <-
    ComputeLayout(
      slot(object, name = "cellgraph"),
      layout_method = layout_method,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      pivots = pivots,
      seed = seed,
      custom_layout_function = custom_layout_function,
      custom_layout_function_args = custom_layout_function_args,
      ...
    )
  if (is.null(slot(object, name = "layout"))) {
    slot(object, name = "layout") <- list()
  }

  # Add layout to CellGraph layout slot
  slot(object, name = "layout")[[layout_name]] <- layout

  return(object)
}

#' @param cl An integer to indicate number of child-processes (integer values are ignored
#' on Windows) for parallel evaluations. See Details on performance
#' in the documentation for \code{pbapply}. The default is NULL,
#' which means that no parallelization is used.
#' @param verbose Print messages
#'
#' @rdname ComputeLayout
#' @method ComputeLayout MPXAssay
#'
#' @examples
#'
#' # Compute layout for a CellGraphAssay
#' cg_assay <- ComputeLayout(seur[["mpxCells"]], layout_method = "wpmds")
#'
#' @export
#'
ComputeLayout.MPXAssay <- function(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
) {
  # Check cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("No 'cellgraphs' loaded. Returning unmodified object.")
    }
    return(object)
  }

  # Only keep loaded graphs
  cellgraphs_loaded <- cellgraphs[loaded_graphs]

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Computing layouts for {length(cellgraphs_loaded)} graphs")
  }

  # Calculate layout for each cell graph object sequentially
  cellgraphs_loaded <- pblapply(cellgraphs_loaded, function(g) {
    g <- ComputeLayout(
      g,
      layout_method = layout_method,
      layout_name = layout_name,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      pivots = pivots,
      seed = seed,
      custom_layout_function = custom_layout_function,
      custom_layout_function_args = custom_layout_function_args,
      ...
    )
    return(g)
  }, cl = cl)

  slot(object, name = "cellgraphs")[names(cellgraphs_loaded)] <- cellgraphs_loaded

  return(object)
}

#' @rdname ComputeLayout
#' @method ComputeLayout CellGraphAssay
#' @docType methods
#' @export
#'
ComputeLayout.CellGraphAssay <- ComputeLayout.MPXAssay

#' @rdname ComputeLayout
#' @method ComputeLayout CellGraphAssay5
#' @docType methods
#' @export
#'
ComputeLayout.CellGraphAssay5 <- ComputeLayout.MPXAssay

#' @rdname ComputeLayout
#' @method ComputeLayout PNAAssay
#'
#' @export
#'
ComputeLayout.PNAAssay <- function(
  object,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
) {
  # Check cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("No 'cellgraphs' loaded. Returning unmodified object.")
    }
    return(object)
  }

  # Only keep loaded graphs
  cellgraphs_loaded <- cellgraphs[loaded_graphs]

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Computing layouts for {length(cellgraphs_loaded)} graphs")
  }

  # Calculate layout for each cell graph object sequentially
  cellgraphs_loaded <- pblapply(cellgraphs_loaded, function(g) {
    g <- ComputeLayout(
      g,
      layout_method = layout_method,
      layout_name = layout_name,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      pivots = pivots,
      seed = seed,
      custom_layout_function = custom_layout_function,
      custom_layout_function_args = custom_layout_function_args,
      ...
    )
    return(g)
  }, cl = cl)

  slot(object, name = "cellgraphs")[names(cellgraphs_loaded)] <- cellgraphs_loaded

  return(object)
}

#' @rdname ComputeLayout
#' @method ComputeLayout PNAAssay5
#' @docType methods
#' @export
#'
ComputeLayout.PNAAssay5 <- ComputeLayout.PNAAssay


#' @param assay Name of assay to compute layouts for
#' @rdname ComputeLayout
#' @method ComputeLayout Seurat
#'
#' @export
#'
ComputeLayout.Seurat <- function(
  object,
  assay = NULL,
  layout_method = c("pmds", "wpmds"),
  layout_name = NULL,
  dim = 3,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  custom_layout_function = NULL,
  custom_layout_function_args = NULL,
  cl = NULL,
  ...
) {
  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  pixel_assay <- object[[assay]]
  assert_class(
    pixel_assay,
    classes = c("CellGraphAssay", "CellGraphAssay5", "PNAAssay", "PNAAssay5")
  )

  # Check cellgraphs
  cellgraphs <- slot(pixel_assay, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("No 'cellgraphs' loaded in 'Seurat' object. Returning unmodified 'Seurat' object.")
    }
    return(object)
  }

  pixel_assay <-
    ComputeLayout(
      pixel_assay,
      layout_method = layout_method,
      layout_name = layout_name,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      pivots = pivots,
      seed = seed,
      verbose = verbose,
      custom_layout_function = custom_layout_function,
      custom_layout_function_args = custom_layout_function_args,
      cl = cl,
      ...
    )

  object[[assay]] <- pixel_assay
  return(object)
}

#' Validate a custom layout function
#'
#' @import rlang
#'
#' @noRd
#'
.validate_custom_layout_function <- function(
  custom_layout_function,
  custom_layout_function_args,
  dim,
  call = caller_env()
) {
  # Make sure that a function is provided
  assert_function(custom_layout_function, call = call)

  # Make sure that the custom_layout_function_args is a list
  if (!is.null(custom_layout_function_args)) {
    if (!inherits(custom_layout_function_args, what = "list")) {
      cli::cli_abort(
        c("x" = "{.var custom_layout_function_args} should be a list with function arguments"),
        call = call
      )
    }
  }
}

#' Validate the results from a custom_layout_function
#'
#' @noRd
#'
.validate_custom_layout_function_results <- function(
  layout,
  dim,
  n_nodes,
  call = caller_env()
) {
  # Validate results
  if (!inherits(layout, what = "matrix")) {
    cli::cli_abort(
      c("x" = "Expected a {.cls matrix} from 'custom_layout_function', but got a {.cls {class(layout)[1]}}"),
      call = call
    )
  }
  if (ncol(layout) != dim) {
    cli::cli_abort(
      c("x" = "Expected {dim} columns in the output from 'custom_layout_function', but got {ncol(layout)}"),
      call = call
    )
  }
  if (nrow(layout) != n_nodes) {
    cli::cli_abort(
      c(
        "x" = "Invalid number of rows returned by 'custom_layout_function'. ",
        " " = "The number of rows in the layout should match the number of ",
        " " = "nodes in the graph ({.val {n_nodes}})"
      ),
      call = call
    )
  }
}


#' Layout Coordinates Utility Functions
#'
#' Utility function used to manipulate layout coordinates.
#' These functions always returns a \code{tbl_df} object.
#'
#' @section Center Layout Coordinates:
#' Centers each axis of the layout coordinates around their means.
#'
#' @param layout A matrix-like object with layout coordinates
#'
#' @name layout coordinates utils
#' @rdname layout_coordinates_utils
#'
#' @return A \code{tbl_df} object with adjusted layout coordinates
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#'
#' # Generate random points that are offset to (100, 100, 100)
#' xyz <- matrix(rnorm(600, mean = 100, sd = 20),
#'   ncol = 3,
#'   dimnames = list(NULL, c("x", "y", "z"))
#' ) %>%
#'   as_tibble()
#'
#' # Visualize random points
#' plotly::plot_ly(
#'   data = xyz, x = ~x, y = ~y, z = ~z,
#'   type = "scatter3d", mode = "markers"
#' )
#'
#' # Center points at (0, 0, 0)
#' xyz_centered <- center_layout_coordinates(xyz)
#' apply(xyz_centered, 2, mean)
#' plotly::plot_ly(
#'   data = xyz_centered, x = ~x, y = ~y,
#'   z = ~z, type = "scatter3d", mode = "markers"
#' )
#'
#' # Normalize points to have a median radius of 1
#' xyz_normalized <- normalize_layout_coordinates(xyz_centered)
#' radii <- sqrt(rowSums(xyz_normalized^2))
#' median(radii)
#' plotly::plot_ly(
#'   data = xyz_normalized, x = ~x, y = ~y,
#'   z = ~z, type = "scatter3d", mode = "markers"
#' )
#'
#' # Project points on unit sphere
#' xyz_projected <- project_layout_coordinates_on_unit_sphere(xyz_normalized)
#' radii <- sqrt(rowSums(xyz_projected^2))
#' all(near(radii, y = rep(1, length(radii)), tol = 1e-12))
#' plotly::plot_ly(
#'   data = xyz_projected, x = ~x, y = ~y,
#'   z = ~z, type = "scatter3d", mode = "markers"
#' )
#'
#' @export
#'
center_layout_coordinates <- function(
  layout
) {
  assert_non_empty_object(layout, classes = c("matrix", "data.frame"))
  if (!(ncol(layout) %in% c(2, 3))) {
    cli::cli_abort(
      c(
        "i" = "The layout must have 2 or 3 columns",
        "x" = "You've provided a layout with {.val {ncol(layout)}} columns"
      )
    )
  }

  # Force rename columns to x, y, z
  colnames(layout) <- c("x", "y", "z")[seq_len(ncol(layout))]
  if (inherits(layout, what = "matrix")) {
    layout <- as_tibble(layout)
  }

  # Center layout coordinates
  layout <- layout %>%
    mutate(across(contains(c("x", "y", "z")), ~ .x - mean(.x)))

  return(layout)
}


#' @section Normalize Layout Coordinates:
#' Centers each axis of the layout coordinates around the mean and
#' adjusts each point coordinate such that their median length
#' (euclidean norm) is 1.
#'
#' @param layout A matrix-like object with layout coordinates
#'
#' @rdname layout_coordinates_utils
#'
#' @export
#'
normalize_layout_coordinates <- function(
  layout
) {
  assert_non_empty_object(layout, classes = c("matrix", "data.frame"))
  if (!(ncol(layout) %in% c(2, 3))) {
    cli::cli_abort(
      c(
        "i" = "The layout must have 2 or 3 columns",
        "x" = "You've provided a layout with {.val {ncol(layout)}} columns"
      )
    )
  }

  layout <- center_layout_coordinates(layout)

  radii <- layout %>%
    mutate(across(contains(c("x", "y", "z")), ~ .x^2)) %>%
    rowSums()
  radii <- sqrt(radii)
  # nolint start
  median_radius <- median(radii)
  layout <- layout %>%
    mutate(across(contains(c("x", "y", "z")), ~ .x / median_radius))
  # nolint end

  return(layout)
}


#' @section Project Layout Coordinates on a Unit Sphere:
#' Centers each axis of the layout coordinates around the mean and
#' adjusts each coordinate such that their lengths (euclidean norm)
#' are 1. This function only accepts layouts with 3 dimensions.
#'
#' @param layout A matrix-like object with layout coordinates
#'
#' @rdname layout_coordinates_utils
#'
#' @export
#'
project_layout_coordinates_on_unit_sphere <- function(
  layout
) {
  assert_non_empty_object(layout, classes = c("matrix", "data.frame"))
  if (ncol(layout) != 3) {
    cli::cli_abort(
      c(
        "i" = "The layout must have 3 columns",
        "x" = "You've provided a layout with {.val {ncol(layout)}} columns"
      )
    )
  }

  layout <- center_layout_coordinates(layout)

  radii <- layout %>%
    mutate(across(contains(c("x", "y", "z")), ~ .x^2)) %>%
    rowSums()
  radii <- sqrt(radii)
  layout <- layout %>%
    mutate(across(contains(c("x", "y", "z")), ~ .x / radii))

  return(layout)
}
