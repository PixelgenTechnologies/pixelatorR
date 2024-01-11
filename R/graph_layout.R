# Declarations used in package check
globalVariables(
  names = c('x', 'y', 'z', 'name'),
  package = 'pixelatorR',
  add = TRUE
)

#' Get a graph layout
#'
#' Calculates a graph layout for a component's edgelist, and outputs a list
#' with the bipartite graph, layout, and antibody counts per A pixel.
#'
#' @param layout_method The method for calculating the graph layout; PMDS, Fruchterman-Reingold (fr), Kamada-Kawai (kk), drl.
#' @param dim An integer specifying the dimensions of the layout (2 or 3)
#' @param normalize_layout Logical specifying whether the coordinate system for the layout should be scaled from -1 to 1.
#' @param k The size of the neighborhood from which to pool counts from in the UPIA antibody count table. 0 is recommended.
#' @param pivots Only used for "pmds" graph layout. See \code{?layout_with_pmds} for details
#' @param project_on_unit_sphere Should the resulting layout be projected onto a unit sphere?
#' @param seed Set seed for reproducibility
#'
#' @import dplyr
#' @importFrom igraph layout_with_fr layout_with_kk layout_with_drl
#'
#' @rdname ComputeLayout
#' @method ComputeLayout tbl_graph
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#'
#' # Load example data
#' seur <- ReadMPX_Seurat(pxl_file, return_cellgraphassay = TRUE, overwrite = TRUE)
#'
#' # Load 1 cellgraph
#' seur <- LoadCellGraphs(seur, cells = colnames(seur)[1],
#'                        load_as = "Anode", force = TRUE)
#'
#' # Get CellGraph
#' cg <- CellGraphs(seur)[[colnames(seur)[1]]]
#'
#' # Get tbl_graph object
#' tbl_graph <- slot(cg, name = "cellgraph")
#'
#' # Compute layout for a tbl_graph
#' layout <- ComputeLayout(tbl_graph, layout_method = "fr")
#' layout %>% head()
#'
#' @export
#'
ComputeLayout.tbl_graph <- function (
  object,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  ...
) {

  # Validate input parameters
  stopifnot("dim should be an integer vector of length 1" = inherits(dim, what = c("numeric", "integer")) & (length(dim) == 1),
            "normalize_layout should be TRUE/FALSE" = is.logical(normalize_layout),
            "k should set to 1 or 2" = inherits(k, what = c("numeric", "integer")) & (length(dim) == 1) & (dim %in% c(2, 3)),
            "'seed' should be a numeric value" = inherits(seed, what = "numeric"))

  # Check and select a layout method
  layout_method <- match.arg(layout_method, choices = c("pmds", "fr", "kk", "drl"))

  layout_function <- switch(
    layout_method,
    "fr" = layout_with_fr,
    "kk" = layout_with_kk,
    "drl" = layout_with_drl,
    "pmds" = {
      # Abort if pmds is selected and graphlayouts isn't installed
      expect_graphlayouts()
      graphlayouts::layout_with_pmds
    })

  # Set seed
  old_seed <- .Random.seed
  set.seed(seed)

  layout <-
    object %>%
    {
      if (layout_method == "pmds") {
        layout_function(., dim = dim, pivots = pivots)
      } else {
        layout_function(., dim = dim, ...)
      }
    } %>%
    as_tibble(.name_repair = function(x) c("x", "y", "z")[1:length(x)]) %>%
    bind_cols(as_tibble(object)) %>%
    select(any_of(c("x", "y", "z")))

  # Project to unit sphere
  if (project_on_unit_sphere || normalize_layout) {

    layout <- layout %>%
      mutate(across(contains(c("x", "y", "z")), ~ .x - mean(x)))

    radii <- layout %>%
      mutate(across(contains(c("x", "y", "z")), ~ .x^2)) %>%
      rowSums() %>%
      sqrt()
    max_radius <- max(radii)
    if (project_on_unit_sphere) {
      layout <- layout %>%
        mutate(across(contains(c("x", "y", "z")), ~ (.x*max_radius)/radii))
    }
    if (normalize_layout) {
      layout <- layout %>%
        mutate(across(contains(c("x", "y", "z")), ~ .x/max_radius))
    }
  }

  # Restore old seed
  .Random.seed <- old_seed

  return(layout)
}


#' @rdname ComputeLayout
#' @method ComputeLayout CellGraph
#'
#' @importFrom future.apply future_lapply
#'
#' @examples
#'
#' # Compute layout for a CellGraph
#' cg <- ComputeLayout(cg, layout_method = "fr")
#'
#' @export
#'
ComputeLayout.CellGraph <- function (
  object,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  ...
) {
  layout <-
    ComputeLayout(
      slot(object, name = "cellgraph"),
      layout_method = layout_method,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      k = k,
      pivots = pivots,
      seed = seed
    )
  if (is.null(slot(object, name = "layout"))) {
    slot(object, name = "layout") <- list()
  }
  slot(object, name = "layout")[[match.arg(layout_method, choices = c("pmds", "fr", "kk", "drl"))]] <- layout
  return(object)
}


#' @param verbose Print messages
#'
#' @rdname ComputeLayout
#' @method ComputeLayout CellGraphAssay
#'
#'
#' @examples
#'
#' # Compute layout for a CellGraphAssay
#' cg_assay <- ComputeLayout(seur[["mpxCells"]], layout_method = "fr")
#'
#' @export
#'
ComputeLayout.CellGraphAssay <- function (
  object,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  ...
) {

  # Check cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity())
      cli_alert_info("No 'cellgraphs' loaded in 'CellGraphAssay'. Returning unmodified 'CellGraphAssay'.")
    return(object)
  }

  # Only keep loaded graphs
  cellgraphs_loaded <- cellgraphs[loaded_graphs]

  if (verbose & check_global_verbosity())
    cli_alert_info("Computing layouts for {length(cellgraphs_loaded)} graphs")

  # Calculate layout for each cell graph object sequentially
  p <- progressr::progressor(along = cellgraphs_loaded)
  cellgraphs_loaded <- future_lapply(cellgraphs_loaded, function(g) {
    g <- ComputeLayout(
      g,
      layout_method = layout_method,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      k = k,
      pivots = pivots,
      seed = seed
    )
    p()
    return(g)
  }, future.seed = seed)

  slot(object, name = "cellgraphs")[names(cellgraphs_loaded)] <- cellgraphs_loaded

  return(object)
}

#' @param assay Name of assay to compute layouts for
#' @rdname ComputeLayout
#' @method ComputeLayout Seurat
#'
#' @examples
#'
#' # Compute layout for a Seurat object
#' seur <- ComputeLayout(seur, layout_method = "fr")
#'
#' @export
#'
ComputeLayout.Seurat <- function (
  object,
  assay = NULL,
  layout_method = c("pmds", "fr", "kk", "drl"),
  dim = 2,
  normalize_layout = FALSE,
  project_on_unit_sphere = FALSE,
  k = 0,
  pivots = 100,
  seed = 123,
  verbose = TRUE,
  ...
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("assay '{assay}' is not a 'CellGraphAssay'"))
  }

  # Check cellgraphs
  cellgraphs <- slot(cg_assay, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity())
      cli_alert_info("No 'cellgraphs' loaded in 'Seurat' object. Returning unmodified 'Seurat' object.")
    return(object)
  }

  cg_assay <-
    ComputeLayout(
      cg_assay,
      layout_method = layout_method,
      dim = dim,
      normalize_layout = normalize_layout,
      project_on_unit_sphere = project_on_unit_sphere,
      k = k,
      pivots = pivots,
      seed = seed,
      verbose = verbose,
      ...
    )

  object[[assay]] <- cg_assay
  return(object)
}

