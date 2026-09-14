#' @include generics.R
NULL

#' @param layout_method Name of a stored layout, typically one computed with
#' \code{\link{ComputeLayout}} or loaded with \code{\link{LoadCellGraphs}}.
#' Default is \code{"wpmds_3d"}. The layout must contain \code{x}, \code{y},
#' and \code{z} coordinates.
#' @param vars Optional character vector of node-level variables to fetch with
#' \code{\link[SeuratObject]{FetchData}} (markers, metadata columns, graph
#' vertex attributes, or reduction embeddings). A variable that is present
#' on some graphs and missing on others is filled with \code{NA}. A variable
#' missing from every graph is omitted, with the same warning as
#' \code{\link[SeuratObject]{FetchData}}.
#' @param add_protein If \code{TRUE}, add a \code{protein} column with the
#' marker label of each node. Labels are read from the protein counts matrix. 
#' Nodes with no count are \code{NA}.
#' @param layer Name of a node matrix layer passed to
#' \code{\link[SeuratObject]{FetchData}}. \code{NULL} (default) uses the same
#' layer selection as \code{FetchData.CellGraph}.
#'
#' @rdname FetchLayoutData
#' @method FetchLayoutData CellGraph
#'
#' @examples
#' library(pixelatorR)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1], add_layouts = TRUE, verbose = FALSE)
#' cg <- CellGraphs(se)[[1]]
#'
#' # Coordinates only
#' layout <- FetchLayoutData(cg)
#'
#' # Include marker counts
#' layout <- FetchLayoutData(cg, vars = "B2M")
#'
#' # Include protein labels
#' layout <- FetchLayoutData(cg, add_protein = TRUE)
#'
#' @export
#'
FetchLayoutData.CellGraph <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  add_protein = FALSE,
  layer = NULL,
  ...
) {
  .fetch_layout_data_graph(
    object,
    layout_method = layout_method,
    vars = vars,
    add_protein = add_protein,
    layer = layer,
    warn_missing = TRUE
  )
}

#' @param cells Component IDs to fetch. If \code{NULL}, all loaded
#' \code{\link{CellGraph}} objects are used. Unloaded graphs raise an error
#' when they are included in \code{cells}.
#'
#' @rdname FetchLayoutData
#' @method FetchLayoutData CellGraphList
#'
#' @examples
#' # Combine layouts from a CellGraphList
#' cgl <- CellGraphs(se)
#' layout <- FetchLayoutData(cgl, cells = colnames(se)[1], vars = "B2M")
#'
#' @export
#'
FetchLayoutData.CellGraphList <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  add_protein = FALSE,
  layer = NULL,
  ...
) {
  cells <- .resolve_fetch_layout_data_cells(object, cells)

  fetched <- dplyr::bind_rows(lapply(cells, function(nm) {
    .fetch_layout_data_graph(
      object[[nm]],
      layout_method = layout_method,
      vars = vars,
      add_protein = add_protein,
      layer = layer,
      warn_missing = FALSE
    ) %>%
      mutate(component = nm, .before = 1)
  }))
  .warn_unfound_fetch_vars(vars, names(fetched))
  fetched
}

#' @rdname FetchLayoutData
#' @method FetchLayoutData PNAAssay
#'
#' @examples
#' # PNAAssay method
#' layout <- FetchLayoutData(se[["PNA"]], cells = colnames(se)[1], vars = "B2M")
#'
#' @export
#'
FetchLayoutData.PNAAssay <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  add_protein = FALSE,
  layer = NULL,
  ...
) {
  FetchLayoutData(
    CellGraphs(object),
    layout_method = layout_method,
    vars = vars,
    cells = cells,
    add_protein = add_protein,
    layer = layer,
    ...
  )
}

#' @rdname FetchLayoutData
#' @method FetchLayoutData PNAAssay5
#' @docType methods
#' @export
#'
FetchLayoutData.PNAAssay5 <- FetchLayoutData.PNAAssay

#' @param assay Name of assay to fetch layouts from
#'
#' @rdname FetchLayoutData
#' @method FetchLayoutData Seurat
#'
#' @examples
#' # Seurat method
#' layout <- FetchLayoutData(se, cells = colnames(se)[1], vars = "B2M")
#'
#' @export
#'
FetchLayoutData.Seurat <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  assay = NULL,
  add_protein = FALSE,
  layer = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  pixel_assay <- object[[assay]]
  assert_pixel_assay(pixel_assay)

  FetchLayoutData(
    CellGraphs(pixel_assay),
    layout_method = layout_method,
    vars = vars,
    cells = cells,
    add_protein = add_protein,
    layer = layer,
    ...
  )
}

#' Fetch a named 3D layout from a CellGraph
#'
#' Looks up \code{layout_method} in the \code{layout} slot and checks that
#' the table has \code{x}, \code{y}, and \code{z} columns.
#'
#' @param object A \code{CellGraph}
#' @param layout_method Name of a stored layout (for example \code{"wpmds_3d"})
#' @param call Environment to report as the error caller
#'
#' @return A data frame of coordinates
#'
#' @keywords internal
#' @noRd
#'
.get_cellgraph_layout <- function(object, layout_method, call = caller_env()) {
  layouts <- slot(object, name = "layout")
  if (is.null(layouts) || length(layouts) == 0) {
    cli::cli_abort(
      c(
        "x" = "This {.cls CellGraph} has no layouts.",
        "i" = "Compute one with {.fn ComputeLayout} or load precomputed layouts with ",
        " " = "{.code LoadCellGraphs(..., add_layouts = TRUE)}."
      ),
      call = call
    )
  }
  if (!layout_method %in% names(layouts)) {
    cli::cli_abort(
      c(
        "x" = "Missing layout {.str {layout_method}}.",
        "i" = "Available layout{?s}: {.val {names(layouts)}}"
      ),
      call = call
    )
  }
  layout <- layouts[[layout_method]]
  if (!inherits(layout, "data.frame")) {
    layout <- as.data.frame(layout, stringsAsFactors = FALSE, check.names = FALSE)
  }
  missing_coords <- setdiff(c("x", "y", "z"), colnames(layout))
  if (length(missing_coords) > 0) {
    cli::cli_abort(
      c(
        "x" = "Layout {.str {layout_method}} must contain {.val x}, {.val y}, and {.val z} coordinates.",
        "i" = "Missing column{?s}: {.val {missing_coords}}"
      ),
      call = call
    )
  }
  layout
}

#' Protein labels from a one-hot encoded node counts matrix
#'
#' Each node has a single protein. The counts matrix is one-hot encoded, so
#' the label is the column name of the non-zero entry in that row. The
#' compressed sparse column index is used so the dense matrix is never
#' materialized.
#'
#' @param object A \code{CellGraph}
#' @param nodes Node names to return labels for, in that order
#' @param call Environment to report as the error caller
#'
#' @return A character vector of protein names, with \code{NA} for nodes
#' that have no count
#'
#' @keywords internal
#' @noRd
#'
.node_protein_labels <- function(object, nodes, call = caller_env()) {
  node_map <- .cg_node_map(object)
  labels <- rep(NA_character_, length(nodes))
  counts <- slot(object, "counts")
  if (is.null(counts) || nrow(counts) == 0L || ncol(counts) == 0L) {
    return(labels)
  }
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }
  if (length(counts@i) == 0L) {
    return(labels)
  }
  proteins <- colnames(counts)
  if (is.null(proteins)) {
    cli::cli_abort(
      c("x" = "The counts matrix has no column names to use as protein labels."),
      call = call
    )
  }
  row_idx <- counts@i + 1L
  col_idx <- rep.int(seq_len(ncol(counts)), diff(counts@p))
  all_labels <- rep(NA_character_, nrow(counts))
  all_labels[row_idx] <- proteins[col_idx]
  all_labels[.row_index(nodes, node_map)]
}

#' Fetch a 3D layout and optional node variables from one CellGraph
#'
#' @param object A \code{CellGraph}
#' @param layout_method Name of a stored layout
#' @param vars Character vector of variable names, or \code{NULL}
#' @param add_protein If \code{TRUE}, add a \code{protein} column
#' @param layer Layer name passed to \code{FetchData}, or \code{NULL}
#' @param warn_missing If \code{TRUE}, warn for requested variables that
#' were not found on this graph
#'
#' @return A tibble of coordinates and fetched variables
#'
#' @keywords internal
#' @noRd
#'
.fetch_layout_data_graph <- function(
  object,
  layout_method,
  vars,
  add_protein,
  layer,
  warn_missing
) {
  assert_single_value(layout_method, type = "string")
  assert_vector(vars, type = "character", n = 1, allow_null = TRUE)
  assert_single_value(add_protein, type = "bool")
  assert_single_value(layer, type = "string", allow_null = TRUE)
  .assert_current_cellgraph(object)

  if (!is.null(vars)) {
    vars <- as.character(vars)
    reserved <- intersect(vars, c("x", "y", "z", "component", "protein"))
    if (length(reserved) > 0) {
      cli::cli_abort(
        c("x" = "{.arg vars} cannot include reserved column name{?s} {.val {reserved}}.")
      )
    }
  }

  layout <- .get_cellgraph_layout(object, layout_method)
  node_names <- .explicit_rownames(layout) %||% .cg_node_map(object)

  fetched <- .fetch_layout_vars(
    object = object,
    vars = vars,
    cells = node_names,
    layer = layer,
    fill_missing = FALSE
  )
  if (isTRUE(warn_missing)) {
    .warn_unfound_fetch_vars(vars, names(fetched))
  }

  coords <- as_tibble(layout[, c("x", "y", "z"), drop = FALSE])
  if (isTRUE(add_protein)) {
    coords$protein <- .node_protein_labels(object, nodes = node_names)
  }
  dplyr::bind_cols(coords, as_tibble(fetched, .name_repair = "minimal"))
}

#' Warn for requested variables that were not found
#'
#' Matches the \code{FetchData.CellGraph} warning when some requested
#' names are absent.
#'
#' @param vars Requested variable names, or \code{NULL}
#' @param found Names present in the result (including reserved layout columns)
#'
#' @return \code{NULL}, invisibly
#'
#' @keywords internal
#' @noRd
#'
.warn_unfound_fetch_vars <- function(vars, found) {
  if (is.null(vars) || length(vars) == 0) {
    return(invisible(NULL))
  }
  vars_missing <- setdiff(vars, found)
  if (length(vars_missing) == 0) {
    return(invisible(NULL))
  }
  m2 <- if (length(vars_missing) > 10) {
    paste0(" (10 out of ", length(vars_missing), " shown)")
  } else {
    ""
  }
  cli::cli_warn(
    "The following requested variables were not found{m2}: {.val {head(vars_missing, 10)}}"
  )
  invisible(NULL)
}

#' Fetch vars for layout rows, filling missing values with NA
#'
#' Calls \code{FetchData} on the \code{CellGraph} and aligns the result
#' to \code{cells}. Variables that are missing stay \code{NA} instead of
#' aborting when \code{fill_missing} is \code{TRUE}. Classed columns such
#' as factors are copied with \code{[[} so types are preserved.
#'
#' @param object A \code{CellGraph}
#' @param vars Character vector of variable names, or \code{NULL}
#' @param cells Node names corresponding to layout rows
#' @param layer Layer name passed to \code{FetchData}, or \code{NULL}
#' @param fill_missing If \code{TRUE}, add a \code{NA} column for each
#' requested variable that was not found. If \code{FALSE}, omit those
#' columns so the caller can warn once after combining graphs.
#' @param call Environment to report as the error caller
#'
#' @return A data frame with rows \code{cells} and columns for the
#' requested variables that were found (and, when \code{fill_missing}
#' is \code{TRUE}, \code{NA} columns for the rest)
#'
#' @keywords internal
#' @noRd
#'
.fetch_layout_vars <- function(
  object,
  vars,
  cells,
  layer,
  fill_missing = TRUE,
  call = caller_env()
) {
  fetched <- data.frame(row.names = cells, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(vars) || length(vars) == 0) {
    return(fetched)
  }

  fetched_data <- tryCatch(
    withCallingHandlers(
      FetchData(
        object,
        vars = vars,
        cells = cells,
        layer = layer,
        clean = FALSE
      ),
      warning = function(w) {
        msg <- conditionMessage(w)
        if (grepl("not found|missing data for vars|not present in this", msg)) {
          tryInvokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("Unknown layer|has no layers", msg)) {
        cli::cli_abort(c("x" = "{msg}"), call = call)
      }
      if (grepl("None of the requested variables|None of the requested nodes", msg)) {
        return(data.frame(row.names = cells, stringsAsFactors = FALSE, check.names = FALSE))
      }
      cli::cli_abort(c("x" = "{msg}"), call = call)
    }
  )

  for (v in vars) {
    if (v %in% colnames(fetched_data) && nrow(fetched_data) > 0) {
      idx <- match(cells, rownames(fetched_data))
      fetched[[v]] <- fetched_data[[v]][idx]
    } else if (isTRUE(fill_missing)) {
      fetched[[v]] <- NA
    }
  }
  keep <- if (isTRUE(fill_missing)) vars else intersect(vars, names(fetched))
  fetched[, keep, drop = FALSE]
}

#' Resolve component IDs for FetchLayoutData
#'
#' Forwards to \code{\link{.resolve_loaded_cellgraph_ids}} with \code{fn = "FetchLayoutData"}.
#'
#' @param object A \code{CellGraphList} (or named list of graphs)
#' @param cells Component IDs to keep, or \code{NULL} for all loaded graphs
#' @param call Environment to report as the error caller
#'
#' @return A character vector of component IDs
#'
#' @keywords internal
#' @noRd
#'
.resolve_fetch_layout_data_cells <- function(object, cells, call = caller_env()) {
  .resolve_loaded_cellgraph_ids(object, cells, fn = "FetchLayoutData", call = call)
}

#' Resolve component IDs and require loaded CellGraph objects
#'
#' When \code{cells} is \code{NULL}, all loaded graphs in the
#' \code{CellGraphList} are used. Requested IDs that are missing or still
#' \code{NULL} abort with a message to run \code{LoadCellGraphs}.
#'
#' @param object A \code{CellGraphList} (or named list of graphs)
#' @param cells Component IDs to keep, or \code{NULL} for all loaded graphs
#' @param fn Name of the calling function, used in error messages
#' @param call Environment to report as the error caller
#'
#' @return A character vector of component IDs
#'
#' @keywords internal
#' @noRd
#'
.resolve_loaded_cellgraph_ids <- function(object, cells, fn, call = caller_env()) {
  available <- names(object)
  loaded <- vapply(object, function(x) inherits(x, "CellGraph"), logical(1))
  loaded_ids <- available[loaded]

  if (is.null(cells)) {
    if (length(loaded_ids) == 0) {
      cli::cli_abort(
        c(
          "x" = "No {.cls CellGraph} objects are loaded.",
          "i" = "Load them with {.fn LoadCellGraphs} before calling {.fn {fn}}."
        ),
        call = call
      )
    }
    return(loaded_ids)
  }

  assert_vector(cells, type = "character", n = 1, call = call)
  assert_unique(cells, call = call)
  assert_x_in_y(cells, available, call = call)

  not_loaded <- cells[!cells %in% loaded_ids]
  if (length(not_loaded) > 0) {
    cli::cli_abort(
      c(
        "x" = "{cli::qty(length(not_loaded))}{.cls CellGraph} object{?s} {?is/are} not
                loaded for {length(not_loaded)} component{?s}: {.val {head(not_loaded, 5)}}.",
        "i" = "Load them with {.fn LoadCellGraphs} before calling {.fn {fn}}."
      ),
      call = call
    )
  }
  cells
}
