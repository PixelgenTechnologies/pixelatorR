#' @include generics.R
NULL

#' @param layout_method Name of a stored layout, typically one computed with
#' \code{\link{ComputeLayout}} or loaded with \code{\link{LoadCellGraphs}}.
#' Default is \code{"wpmds_3d"}. The layout must contain \code{x}, \code{y},
#' and \code{z} coordinates.
#' @param vars Optional character vector of node-level variables to fetch with
#' \code{\link[SeuratObject]{FetchData}} (markers, metadata columns, graph
#' vertex attributes, or reduction embeddings). Missing values are filled with
#' \code{NA} rather than raising an error.
#' @param layer Name of a node matrix layer passed to
#' \code{\link[SeuratObject]{FetchData}}. \code{NULL} (default) uses the same
#' layer selection as \code{FetchData.CellGraph}.
#'
#' @rdname ExtractLayout
#' @method ExtractLayout CellGraph
#'
#' @examples
#' library(pixelatorR)
#'
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#' se <- LoadCellGraphs(se, cells = colnames(se)[1], add_layouts = TRUE, verbose = FALSE)
#' cg <- CellGraphs(se)[[1]]
#'
#' # Coordinates only
#' layout <- ExtractLayout(cg)
#'
#' # Include marker counts
#' layout <- ExtractLayout(cg, vars = "B2M")
#'
#' @export
#'
ExtractLayout.CellGraph <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  layer = NULL,
  ...
) {
  object <- .upgrade_cellgraph(object)
  assert_single_value(layout_method, type = "string")
  assert_vector(vars, type = "character", n = 1, allow_null = TRUE)
  assert_single_value(layer, type = "string", allow_null = TRUE)

  if (!is.null(vars)) {
    vars <- as.character(vars)
    reserved <- intersect(vars, c("x", "y", "z", "component"))
    if (length(reserved) > 0) {
      cli::cli_abort(
        c("x" = "{.arg vars} cannot include reserved column name{?s} {.val {reserved}}.")
      )
    }
  }

  layout <- .get_cellgraph_layout(object, layout_method)
  node_names <- .explicit_rownames(layout) %||% Cells(object)

  fetched <- .fetch_layout_vars(
    object = object,
    vars = vars,
    cells = node_names,
    layer = layer
  )

  coords <- as_tibble(layout[, c("x", "y", "z"), drop = FALSE])
  dplyr::bind_cols(coords, as_tibble(fetched, .name_repair = "minimal"))
}

#' @param cells Component IDs to extract. If \code{NULL}, all loaded
#' \code{\link{CellGraph}} objects are used. Unloaded graphs raise an error
#' when they are included in \code{cells}.
#'
#' @rdname ExtractLayout
#' @method ExtractLayout CellGraphList
#'
#' @examples
#' # Combine layouts from a CellGraphList
#' cgl <- CellGraphs(se)
#' layout <- ExtractLayout(cgl, cells = colnames(se)[1], vars = "B2M")
#'
#' @export
#'
ExtractLayout.CellGraphList <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  layer = NULL,
  ...
) {
  cells <- .resolve_extract_layout_cells(object, cells)

  dplyr::bind_rows(lapply(cells, function(nm) {
    ExtractLayout(
      object[[nm]],
      layout_method = layout_method,
      vars = vars,
      layer = layer,
      ...
    ) %>%
      mutate(component = nm, .before = 1)
  }))
}

#' @rdname ExtractLayout
#' @method ExtractLayout PNAAssay
#'
#' @examples
#' # PNAAssay method
#' layout <- ExtractLayout(se[["PNA"]], cells = colnames(se)[1], vars = "B2M")
#'
#' @export
#'
ExtractLayout.PNAAssay <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  layer = NULL,
  ...
) {
  ExtractLayout(
    CellGraphs(object),
    layout_method = layout_method,
    vars = vars,
    cells = cells,
    layer = layer,
    ...
  )
}

#' @rdname ExtractLayout
#' @method ExtractLayout PNAAssay5
#' @docType methods
#' @export
#'
ExtractLayout.PNAAssay5 <- ExtractLayout.PNAAssay

#' @param assay Name of assay to extract layouts from
#'
#' @rdname ExtractLayout
#' @method ExtractLayout Seurat
#'
#' @examples
#' # Seurat method
#' layout <- ExtractLayout(se, cells = colnames(se)[1], vars = "B2M")
#'
#' @export
#'
ExtractLayout.Seurat <- function(
  object,
  layout_method = "wpmds_3d",
  vars = NULL,
  cells = NULL,
  assay = NULL,
  layer = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object)
  pixel_assay <- object[[assay]]
  assert_pixel_assay(pixel_assay)

  ExtractLayout(
    CellGraphs(pixel_assay),
    layout_method = layout_method,
    vars = vars,
    cells = cells,
    layer = layer,
    ...
  )
}

#' Fetch a named 3D layout from a CellGraph
#'
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

#' Fetch vars for layout rows, filling missing values with NA
#'
#' @noRd
#'
.fetch_layout_vars <- function(object, vars, cells, layer, call = caller_env()) {
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
    fetched[[v]] <- NA
  }
  if (ncol(fetched_data) > 0 && nrow(fetched_data) > 0) {
    common_cells <- intersect(rownames(fetched_data), cells)
    common_vars <- intersect(colnames(fetched_data), vars)
    if (length(common_cells) > 0 && length(common_vars) > 0) {
      fetched[common_cells, common_vars] <- fetched_data[common_cells, common_vars, drop = FALSE]
    }
  }
  fetched[, vars, drop = FALSE]
}

#' Resolve component IDs and require loaded CellGraph objects
#'
#' @noRd
#'
.resolve_extract_layout_cells <- function(object, cells, call = caller_env()) {
  available <- names(object)
  loaded <- vapply(object, function(x) inherits(x, "CellGraph"), logical(1))
  loaded_ids <- available[loaded]

  if (is.null(cells)) {
    if (length(loaded_ids) == 0) {
      cli::cli_abort(
        c(
          "x" = "No {.cls CellGraph} objects are loaded.",
          "i" = "Load them with {.fn LoadCellGraphs} before calling {.fn ExtractLayout}."
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
        "x" = "{.cls CellGraph} object{?s} {?is/are} not loaded for {length(not_loaded)} component{?s}: {.val {head(not_loaded, 5)}}.",
        "i" = "Load them with {.fn LoadCellGraphs} before calling {.fn ExtractLayout}."
      ),
      call = call
    )
  }
  cells
}
