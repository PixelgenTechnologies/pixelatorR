#' @include generics.R local_proximity_scores.R CellGraphList.R
NULL

#' Compute local proximity scores
#'
#' Computes local proximity scores for nodes in PNA cell graphs using
#' \code{\link{local_proximity}} and stores the output on each
#' \code{\link{CellGraph}}.
#'
#' The default \code{mode} is \code{"self-clustering"}, which returns a
#' node-by-marker matrix that is stored as a layer. Other modes return a
#' single score per node, which is stored in \code{meta.data}.
#'
#' @param markers A character vector specifying the markers to use. If
#' \code{NULL}, all markers in the count matrix of each \code{CellGraph}
#' are used. Methods that iterate over multiple \code{CellGraph} objects
#' keep the intersection with available markers. If none of the requested
#' markers are present in a graph, a warning is emitted and that graph is
#' left unmodified.
#' @param method A character string specifying the method to use for
#' computing the local proximity score. Options are \code{"analytical"}
#' or \code{"permutation"}.
#' @param mode A character string specifying the mode of computation.
#' See \code{\link{local_proximity}} for details. Default is
#' \code{"self-clustering"}.
#' @param iterations An integer specifying the number of iterations to run
#' when \code{method = "permutation"}.
#' @param k An integer specifying the neighborhood size to consider.
#' @param A_k An optional pre-computed expanded adjacency matrix. Only
#' used by the \code{CellGraph} method.
#' @param seed An integer for random seed setting.
#' @param name Name of the layer (when a matrix is returned) or metadata
#' column (when a vector is returned) used to store scores. Default is
#' \code{"lps"}.
#'
#' @rdname ComputeLPS
#' @method ComputeLPS CellGraph
#'
#' @examples
#' library(pixelatorR)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE)
#' cg <- CellGraphs(se)[[1]]
#'
#' # Matrix result is stored as a layer
#' cg <- ComputeLPS(cg, markers = "B2M")
#' Layers(cg)
#'
#' # Vector result is stored in node metadata
#' cg <- ComputeLPS(cg, markers = "B2M", mode = "all", name = "lps_b2m")
#' head(CellGraphData(cg, slot = "meta.data"))
#'
#' @export
#'
ComputeLPS.CellGraph <- function(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  A_k = NULL,
  seed = 123,
  name = "lps",
  ...
) {
  object <- .upgrade_cellgraph(object)
  method <- match.arg(method, choices = c("analytical", "permutation"))
  mode <- match.arg(mode, choices = c("self-clustering", "all", "any"))
  assert_single_value(name, type = "string")

  counts <- slot(object, "counts")
  if (is.null(counts)) {
    cli::cli_abort(
      c("x" = "{.cls CellGraph} has no counts. Cannot compute local proximity scores.")
    )
  }
  if (is.null(markers)) {
    markers <- colnames(counts)
  }

  scores <- local_proximity(
    object = object,
    markers = markers,
    method = method,
    mode = mode,
    iterations = iterations,
    k = k,
    A_k = A_k,
    seed = seed,
    ...
  )

  if (is.null(dim(scores))) {
    object <- AddMetaData(object, metadata = scores, col.name = name)
  } else {
    if (is.null(colnames(scores)) && !is.null(markers) && ncol(scores) == length(markers)) {
      colnames(scores) <- markers
    }
    LayerData(object, layer = name) <- scores
  }

  return(object)
}

#' @param cl An integer to indicate number of child-processes (integer values
#' are ignored on Windows) for parallel evaluations. See Details on performance
#' in the documentation for \code{pbapply}. The default is \code{NULL},
#' which means that no parallelization is used.
#' @param verbose Print messages
#'
#' @rdname ComputeLPS
#' @method ComputeLPS CellGraphList
#'
#' @export
#'
ComputeLPS.CellGraphList <- function(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
) {
  method <- match.arg(method, choices = c("analytical", "permutation"))
  mode <- match.arg(mode, choices = c("self-clustering", "all", "any"))
  cellgraphs <- slot(object, "cellgraphs")

  if (length(cellgraphs) == 0) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("No CellGraph objects in {.cls CellGraphList}. Returning unmodified object.")
    }
    return(object)
  }

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Computing local proximity scores for {length(cellgraphs)} graph{?s}")
  }

  nms <- names(cellgraphs)
  cellgraphs <- pblapply(nms, function(nm) {
    .compute_lps_cellgraph(
      object = cellgraphs[[nm]],
      markers = markers,
      method = method,
      mode = mode,
      iterations = iterations,
      k = k,
      seed = seed,
      name = name,
      graph_id = nm,
      ...
    )
  }, cl = cl)
  names(cellgraphs) <- nms

  slot(object, "cellgraphs") <- cellgraphs
  return(object)
}

#' @rdname ComputeLPS
#' @method ComputeLPS PNAAssay
#'
#' @examples
#' # Compute LPS for loaded cell graphs in a PNAAssay
#' pna_assay <- ComputeLPS(se[["PNA"]], markers = "B2M")
#'
#' @export
#'
ComputeLPS.PNAAssay <- function(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
) {
  method <- match.arg(method, choices = c("analytical", "permutation"))
  mode <- match.arg(mode, choices = c("self-clustering", "all", "any"))
  cellgraphs <- slot(object, name = "cellgraphs")
  loaded_graphs <- !sapply(cellgraphs, is.null)

  if (sum(loaded_graphs) == 0) {
    if (verbose && check_global_verbosity()) {
      cli_alert_info("No 'cellgraphs' loaded. Returning unmodified object.")
    }
    return(object)
  }

  cellgraphs_loaded <- CreateCellGraphList(cellgraphs[loaded_graphs])
  cellgraphs_loaded <- ComputeLPS(
    cellgraphs_loaded,
    markers = markers,
    method = method,
    mode = mode,
    iterations = iterations,
    k = k,
    seed = seed,
    name = name,
    verbose = verbose,
    cl = cl,
    ...
  )
  slot(object, name = "cellgraphs")[names(cellgraphs_loaded)] <- as.list(cellgraphs_loaded)

  return(object)
}

#' @rdname ComputeLPS
#' @method ComputeLPS PNAAssay5
#' @docType methods
#' @export
#'
ComputeLPS.PNAAssay5 <- ComputeLPS.PNAAssay

#' @param assay Name of assay to compute local proximity scores for
#'
#' @rdname ComputeLPS
#' @method ComputeLPS Seurat
#'
#' @examples
#' # Seurat method (only loaded cell graphs are processed)
#' seur <- ComputeLPS(se, markers = "B2M")
#'
#' @export
#'
ComputeLPS.Seurat <- function(
  object,
  assay = NULL,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  pixel_assay <- object[[assay]]
  assert_pna_assay(pixel_assay)

  pixel_assay <- ComputeLPS(
    pixel_assay,
    markers = markers,
    method = method,
    mode = mode,
    iterations = iterations,
    k = k,
    seed = seed,
    name = name,
    verbose = verbose,
    cl = cl,
    ...
  )

  object[[assay]] <- pixel_assay
  return(object)
}

#' Compute LPS for one CellGraph, intersecting requested markers
#'
#' Used by methods that iterate over many graphs. Missing markers are
#' dropped. If none remain, the graph is returned unmodified with a warning.
#'
#' @noRd
#'
.compute_lps_cellgraph <- function(
  object,
  markers = NULL,
  method,
  mode,
  iterations,
  k,
  seed,
  name,
  graph_id = NULL,
  ...
) {
  if (!is.null(markers)) {
    counts <- slot(object, "counts")
    available <- if (is.null(counts)) character(0) else colnames(counts)
    keep <- intersect(markers, available)
    if (length(keep) == 0) {
      graph_label <- graph_id %||% "CellGraph"
      cli::cli_warn(
        c(
          "!" = "None of the requested markers are present in {.val {graph_label}}.",
          "i" = "Returning the {.cls CellGraph} unmodified."
        )
      )
      return(object)
    }
    markers <- keep
  }

  ComputeLPS(
    object,
    markers = markers,
    method = method,
    mode = mode,
    iterations = iterations,
    k = k,
    seed = seed,
    name = name,
    ...
  )
}
