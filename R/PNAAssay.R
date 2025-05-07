#' @include generics.R lazy_load_tables.R
#' @importClassesFrom Matrix dgCMatrix
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The PNAAssay class
#'
#' The PNAAssay object is an extended \code{\link[SeuratObject]{Assay}}
#' for storing PNA single-cell data.
#'
#' Compared to the \code{\link[SeuratObject]{Assay}} class, the PNAAssay
#' class has three additional slots:
#' - \code{cellgraphs}: A named list of \code{\link{CellGraph}} objects
#' - \code{proximity}: A \code{tbl_df} with proximity scores
#' - \code{fs_map}: A \code{tbl_df} with information on source PXL file
#'
#' @slot cellgraphs A named list of \code{\link{CellGraph}} objects
#' @slot proximity A \code{tbl_df} with proximity scores
#' @slot fs_map A \code{tbl_df} with information on source PXL file
#' paths, sample IDs, and component IDs
#'
#' @name PNAAssay-class
#' @rdname PNAAssay-class
#' @importClassesFrom SeuratObject Assay
#' @exportClass PNAAssay
#' @concept assay
PNAAssay <- setClass(
  Class = "PNAAssay",
  contains = "Assay",
  slots = list(
    cellgraphs = "list",
    proximity = "data.frame",
    fs_map = "data.frame"
  )
)

#' The PNAAssay5 class
#'
#' The PNAAssay5 object is an extended \code{\link[SeuratObject]{Assay5}}
#' for storing PNA single-cell data.
#'
#' Compared to the \code{\link[SeuratObject]{Assay5}} class, the PNAAssay5
#' class has three additional slots:
#' - \code{cellgraphs}: A named list of \code{\link{CellGraph}} objects
#' - \code{proximity}: A \code{tbl_df} with proximity scores
#' - \code{fs_map}: A \code{tbl_df} with information on source PXL file
#'
#' @slot cellgraphs A named list of \code{\link{CellGraph}} objects
#' @slot proximity A \code{tbl_df} with proximity scores
#' @slot fs_map A \code{tbl_df} with information on source PXL file
#' paths, sample IDs, and component IDs
#'
#' @name PNAAssay5-class
#' @rdname PNAAssay5-class
#' @importClassesFrom SeuratObject Assay
#' @exportClass PNAAssay5
#' @concept assay
PNAAssay5 <- setClass(
  Class = "PNAAssay5",
  contains = "Assay5",
  slots = list(
    cellgraphs = "list",
    proximity = "data.frame",
    fs_map = "data.frame"
  )
)



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create a PNAAssay object
#'
#' Create a \code{\link{PNAAssay}} object from a count matrix and a list
#' of \code{\link{CellGraph}} objects. The expected format of the input
#' matrix is features x cells. Optionally, a \code{tbl_df} with proximity
#' scores and a \code{tbl_df} with information on source PXL file paths
#' can be provided.
#'
#' @param counts Unnormalized data (raw counts)
#' @param cellgraphs A named list with \code{\link{CellGraph}} objects
#' @param proximity A \code{tbl_df} with proximity scores
#' @param fs_map A \code{tbl_df} with information on source PXL file
#' paths, sample IDs, and component IDs
#' @param verbose Print messages
#' @param ... Additional arguments passed to \code{\link[SeuratObject]{CreateAssayObject}}
#'
#' @concept assay
#'
#' @return A \code{PNAAssay} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' counts <- ReadPNA_counts(pxl_file)
#' pna_assay <- CreatePNAAssay(
#'   counts = counts,
#'   cellgraphs = rep(list(NULL), ncol(counts)) %>%
#'     setNames(colnames(counts))
#' )
#' pna_assay
#'
#' @export
#'
CreatePNAAssay <- function(
  counts,
  cellgraphs,
  proximity = NULL,
  fs_map = NULL,
  verbose = FALSE,
  ...
) {
  # Check input parameters
  assert_class(counts, c("matrix", "dgCMatrix"))
  assert_class(cellgraphs, "list")
  assert_vectors_match(names(cellgraphs), colnames(counts))

  cellgraphs <- cellgraphs[colnames(counts)]

  # Validate proximity scores
  proximity <-
    .validate_proximity(
      proximity,
      cell_ids = colnames(counts),
      markers = rownames(counts),
      verbose = verbose
    )

  # Create the Seurat assay object
  counts <- as(counts, "dgCMatrix")
  seurat_assay <- CreateAssayObject(
    counts = counts,
    ...
  )

  # Convert Seurat Assay to a CellGraphAssay
  pna_assay <- as.PNAAssay(
    x = seurat_assay,
    cellgraphs = cellgraphs,
    proximity = proximity,
    fs_map = fs_map
  )
  return(pna_assay)
}


#' Create a PNAAssay5 object
#'
#' Create a \code{\link{PNAAssay5}} object from a count matrix and a list
#' of \code{\link{CellGraph}} objects. The expected format of the input
#' matrix is features x cells. Optionally, a \code{tbl_df} with proximity
#' scores and a \code{tbl_df} with information on source PXL file paths
#' can be provided.
#'
#' @param counts Unnormalized data (raw counts)
#' @param cellgraphs A named list with \code{\link{CellGraph}} objects
#' @param proximity A \code{tbl_df} with proximity scores
#' @param fs_map A \code{tbl_df} with information on source PXL file
#' paths, sample IDs, and component IDs
#' @param verbose Print messages
#' @param ... Additional arguments passed to \code{\link[SeuratObject]{CreateAssay5Object}}
#'
#' @concept assay
#'
#' @return A \code{PNAAssay5} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' counts <- ReadPNA_counts(pxl_file)
#' pna_assay5 <- CreatePNAAssay5(
#'   counts = counts,
#'   cellgraphs = rep(list(NULL), ncol(counts)) %>%
#'     setNames(colnames(counts))
#' )
#' pna_assay5
#'
#' @export
#'
CreatePNAAssay5 <- function(
  counts,
  cellgraphs,
  proximity = NULL,
  fs_map = NULL,
  verbose = FALSE,
  ...
) {
  # Check input parameters
  assert_class(counts, c("matrix", "dgCMatrix"))
  assert_class(cellgraphs, "list")
  assert_vectors_match(names(cellgraphs), colnames(counts))

  cellgraphs <- cellgraphs[colnames(counts)]

  # Validate proximity scores
  proximity <-
    .validate_proximity(
      proximity,
      cell_ids = colnames(counts),
      markers = rownames(counts),
      verbose = verbose
    )

  # Create the Seurat assay5 object
  counts <- as(counts, "dgCMatrix")
  seurat_assay <- CreateAssay5Object(
    counts = counts,
    ...
  )

  # Convert Seurat Assay5 to a CellGraphAssay5
  pna_assay5 <- as.PNAAssay5(
    x = seurat_assay,
    cellgraphs = cellgraphs,
    proximity = proximity,
    fs_map = fs_map
  )
  return(pna_assay5)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get/set methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get and set cellgraphs in a PNAAssay or PNAAssay5 object
#'
#' @param object An object with cellgraphs
#' @param value A named list with \code{\link{CellGraph}} objects to
#' replace the current cellgraphs
#' @param ... Additional arguments
#'
#' @rdname CellGraphs
#' @method CellGraphs PNAAssay
#' @export
#' @concept assay
#' @concept cellgraphs
#'
#' @examples
#' library(pixelatorR)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' CellGraphs(seur_obj[["PNA"]])
#'
CellGraphs.PNAAssay <- function(
  object,
  ...
) {
  return(slot(object, name = "cellgraphs"))
}

#' @rdname CellGraphs
#' @method CellGraphs PNAAssay5
#' @export
#' @concept assay
#' @concept cellgraphs
#'
CellGraphs.PNAAssay5 <- CellGraphs.PNAAssay


#' Internal function to replace cellgraphs in a PNAAssay object
#' @noRd
.replace_cellgraphs <- function(
  object,
  value,
  call = caller_env()
) {
  # Clean cellgraphs slot if value = NULL
  if (is.null(x = value)) {
    slot(object = object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% set_names(nm = colnames(object))
    return(object)
  }

  # Validate list
  if (inherits(x = value, what = "list")) {
    assert_vectors_match(names(value), names(CellGraphs(object)))
    for (i in seq_along(value)) {
      if (!inherits(x = value[[i]], what = c("CellGraph", "NULL"))) {
        cli::cli_abort(
          c(
            "i" = "All elements of {.var value} must be of class {.cls CellGraph} or {.cls NULL}",
            "x" = "Element {i} of {.var value} is a {.cls {class(value[[i]])}}"
          ),
          call = call
        )
      }
    }
    slot(object = object, name = "cellgraphs") <- value
  } else {
    cli::cli_abort(
      c(
        "i" = "{.var value} must be a {.cls list}",
        "x" = "You've provided a {.cls {class(value)}}"
      ),
      call = call
    )
  }
  return(object)
}

#' @export
#' @method CellGraphs<- PNAAssay
#' @rdname CellGraphs
#' @concept assay
#' @concept cellgraphs
#'
#' @examples
#' # Set cellgraphs in a PNAAssay object
#' CellGraphs(seur_obj[["PNA"]]) <- CellGraphs(seur_obj[["PNA"]])
#'
"CellGraphs<-.PNAAssay" <- function(
  object,
  ...,
  value
) {
  .replace_cellgraphs(object, value)
}

#' @export
#' @method CellGraphs<- PNAAssay5
#' @rdname CellGraphs
#' @concept assay
#' @concept cellgraphs
#'
"CellGraphs<-.PNAAssay5" <- function(
  object,
  ...,
  value
) {
  .replace_cellgraphs(object, value)
}


#' @param object A \code{PNAAssay} or a \code{PNAAssay5}
#' @param new.names A character vector with new cell IDs. The length of the vector
#' must be equal to the number of cells in the object and the names must be unique.
#' @param ... Additional arguments (not used)
#'
#' @describeIn PNAAssay-methods Rename cell IDs of a \code{PNAAssay} or
#' \code{PNAAssay5} object
#' @method RenameCells PNAAssay
#' @concept assay
#' @docType methods
#'
#' @examples
#' library(pixelatorR)
#' library(SeuratObject)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' pna_assay <- seur_obj[["PNA"]]
#' pna_assay <- RenameCells(pna_assay, new.names = paste0(colnames(pna_assay), "-new"))
#' colnames(pna_assay)
#'
#' @export
RenameCells.PNAAssay <- function(
  object,
  new.names = NULL,
  ...
) {
  assert_vector(new.names, type = "character")
  assert_vectors_x_y_length_equal(new.names, colnames(object))
  assert_unique(new.names)

  if (is.null(names(new.names))) {
    new.names <- set_names(new.names, nm = colnames(object))
  }

  orig.names <- colnames(object)

  assay <- as(object, Class = ifelse(is(object, "PNAAssay"), "Assay", "Assay5"))
  assay <- RenameCells(assay, new.names = new.names, ...)

  # Recreate the CellGraphAssay object
  as_assay_function <- ifelse(is(object, "PNAAssay"), as.PNAAssay, as.PNAAssay5)
  pna_assay_renamed <- as_assay_function(assay)

  # Rename cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  cellgraphs <- set_names(cellgraphs, nm = new.names)
  slot(pna_assay_renamed, name = "cellgraphs") <- cellgraphs

  # Handle proximity slot
  name_conversion <- tibble(new = new.names, component = orig.names)
  proximity <- slot(object, name = "proximity")
  if (length(proximity) > 0) {
    proximity <- proximity %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(pna_assay_renamed, name = "proximity") <- proximity

  # Handle fs_map
  fs_map <- slot(object, name = "fs_map")
  if (nrow(fs_map) > 0) {
    fs_map$id_map <- fs_map %>%
      pull(id_map) %>%
      lapply(function(x) {
        x$current_id <- new.names[x$current_id]
        return(x)
      })
  }
  slot(pna_assay_renamed, name = "fs_map") <- fs_map

  return(pna_assay_renamed)
}

#' @inheritParams RenameCells.PNAAssay
#' @describeIn PNAAssay5-methods Rename cell IDs of a \code{PNAAssay5} object
#' @concept assay
#' @method RenameCells PNAAssay5
#' @docType methods
#' @export
#'
RenameCells.PNAAssay5 <- RenameCells.PNAAssay


#' @param cellgraphs A list of \code{\link{CellGraph}} objects
#' @param proximity A \code{tbl_df} with proximity scores
#' @param fs_map A \code{tbl_df} with information on source PXL file
#' paths, sample IDs, and component IDs
#'
#' @rdname as.PNAAssay
#' @method as.PNAAssay Assay
#' @concept assay
#'
#' @return A \code{PNAAssay} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#' library(SeuratObject)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' counts <- ReadPNA_counts(pxl_file)
#' assay <- CreateAssayObject(
#'   counts = counts
#' )
#' pna_assay <- as.PNAAssay(
#'   assay,
#'   cellgraphs = rep(list(NULL), ncol(counts)) %>%
#'     setNames(colnames(counts))
#' )
#' pna_assay
#'
#' @export
#'
as.PNAAssay.Assay <- function(
  x,
  cellgraphs = NULL,
  proximity = NULL,
  fs_map = NULL,
  ...
) {
  # Check cellgraphs
  if (!is.null(cellgraphs)) {
    assert_non_empty_object(cellgraphs, classes = "list")
    assert_singles_match(length(cellgraphs), ncol(x))
    for (i in seq_along(cellgraphs)) {
      if (!inherits(x = cellgraphs[[i]], what = c("CellGraph", "NULL"))) {
        cli::cli_abort(
          c(
            "i" = "All elements of {.var cellgraphs} must be of class {.cls CellGraph} or {.cls NULL}",
            "x" = "Element {i} of {.var cellgraphs} is a {.cls {class(cellgraphs[[i]])}}"
          )
        )
      }
    }
  } else {
    cellgraphs <- rep(list(NULL), ncol(x)) %>% set_names(nm = colnames(x))
  }

  # Check fs_map
  if (!is.null(fs_map)) {
    .validate_fs_map_pna(fs_map, data_type = "PNA", allow_empty_df = TRUE)
  }

  # Validate proximity scores
  proximity <- .validate_proximity(proximity, cell_ids = colnames(x), markers = rownames(x))

  new_assay <- as(object = x, Class = "PNAAssay")

  # Add slots
  slot(new_assay, name = "cellgraphs") <- cellgraphs
  slot(new_assay, name = "proximity") <- proximity

  # Add fs_map
  if (!is.null(fs_map)) {
    slot(new_assay, name = "fs_map") <- fs_map
  }

  return(new_assay)
}

setAs(
  from = "Assay",
  to = "PNAAssay",
  def = function(from) {
    object.list <- sapply(
      X = slotNames(x = from),
      FUN = slot,
      object = from,
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    object.list <- c(
      list(
        "Class" = "PNAAssay"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)


#' @param cellgraphs A list of \code{\link{CellGraph}} objects
#' @param proximity A \code{tbl_df} with proximity scores
#' @param fs_map A \code{tbl_df} with information on source PXL file
#' paths, sample IDs, and component IDs
#'
#' @rdname as.PNAAssay5
#' @method as.PNAAssay5 Assay5
#' @concept assay
#'
#' @return A \code{PNAAssay5} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#' library(SeuratObject)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' counts <- ReadPNA_counts(pxl_file)
#' assay <- CreateAssay5Object(
#'   counts = counts
#' )
#' pna_assay5 <- as.PNAAssay5(
#'   assay,
#'   cellgraphs = rep(list(NULL), ncol(counts)) %>%
#'     setNames(colnames(counts))
#' )
#' pna_assay5
#'
#' @export
#'
as.PNAAssay5.Assay5 <- function(
  x,
  cellgraphs = NULL,
  proximity = NULL,
  fs_map = NULL,
  ...
) {
  # Check cellgraphs
  if (!is.null(cellgraphs)) {
    assert_non_empty_object(cellgraphs, classes = "list")
    assert_singles_match(length(cellgraphs), ncol(x))
    for (i in seq_along(cellgraphs)) {
      if (!inherits(x = cellgraphs[[i]], what = c("CellGraph", "NULL"))) {
        cli::cli_abort(
          c(
            "i" = "All elements of {.var cellgraphs} must be of class {.cls CellGraph} or {.cls NULL}",
            "x" = "Element {i} of {.var cellgraphs} is a {.cls {class(cellgraphs[[i]])}}"
          )
        )
      }
    }
  } else {
    cellgraphs <- rep(list(NULL), ncol(x)) %>% set_names(nm = colnames(x))
  }

  # Check fs_map
  if (!is.null(fs_map)) {
    .validate_fs_map_pna(fs_map, data_type = "PNA", allow_empty_df = TRUE)
  }

  # Validate proximity scores
  proximity <- .validate_proximity(proximity, cell_ids = colnames(x), markers = rownames(x))

  new_assay <- as(object = x, Class = "PNAAssay5")

  # Add slots
  slot(new_assay, name = "cellgraphs") <- cellgraphs
  slot(new_assay, name = "proximity") <- proximity

  # Add fs_map
  if (!is.null(fs_map)) {
    slot(new_assay, name = "fs_map") <- fs_map
  }

  return(new_assay)
}

setAs(
  from = "Assay5",
  to = "PNAAssay5",
  def = function(from) {
    object.list <- sapply(
      X = slotNames(x = from),
      FUN = slot,
      object = from,
      simplify = FALSE,
      USE.NAMES = TRUE
    )
    object.list <- c(
      list(
        "Class" = "PNAAssay5"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)


#' @param add_marker_counts A logical indicating whether to add marker
#' counts to the output ("count_1" and "count_2")
#' @param add_marker_counts A logical indicating whether to add marker
#' count proportions to the output ("p1" and "p2")
#' @param lazy A logical indicating whether to lazy load the proximity scores
#' from the PXL files
#' @param calc_log2ratio A logical indicating whether to calculate the
#' log2 ratio proximity score
#'
#' @method ProximityScores PNAAssay
#' @rdname ProximityScores
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' ProximityScores(seur_obj[["PNA"]])
#'
#' @export
#'
ProximityScores.PNAAssay <- function(
  object,
  add_marker_counts = FALSE,
  add_marker_proportions = FALSE,
  lazy = FALSE,
  calc_log2ratio = TRUE,
  ...
) {
  assert_single_value(add_marker_counts, type = "bool")
  assert_single_value(add_marker_proportions, type = "bool")
  assert_single_value(calc_log2ratio, type = "bool")
  assert_single_value(lazy, type = "bool")

  if (lazy) {
    # Load proximity scores from the PXL files
    fs_map <- FSMap(object)
    proximity_scores <- .lazy_load_table(fs_map, "proximity", calc_log2ratio, rownames(object))
  } else {
    proximity_scores <- slot(object, name = "proximity")

    if (is.null(proximity_scores) && !lazy) {
      cli::cli_abort(
        c("x" = "proximity scores are missing from the {.cls {class(object)}} object")
      )
    }
  }

  if (add_marker_proportions & !add_marker_counts) {
    cli::cli_alert_warning(
      "Setting {.var add_marker_counts = TRUE} which is required when {.var add_marker_proportions = TRUE}."
    )
    add_marker_counts <- TRUE
  }

  # Add marker counts
  if (add_marker_counts) {
    if (length(proximity_scores) == 0) {
      cli::cli_abort(
        c(
          "i" = "Cannot add marker counts to an empty table",
          "x" = "The proximity score table is empty"
        )
      )
    }
    all_counts <- FetchData(object, vars = rownames(object), layer = "counts") %>%
      rownames_to_column("component") %>%
      pivot_longer(where(is.numeric), names_to = "marker", values_to = "count") %>%
      mutate(count = as.integer(count))

    if (lazy) {
      # Copy to federated database to enable join operation
      copy_to(proximity_scores$src$con, all_counts)
      all_counts <- tbl(proximity_scores$src$con, "all_counts")
    }

    proximity_scores <- proximity_scores %>%
      left_join(all_counts, by = c("component", "marker_1" = "marker")) %>%
      rename(count_1 = count) %>%
      left_join(all_counts, by = c("component", "marker_2" = "marker")) %>%
      rename(count_2 = count)

    if (add_marker_proportions) {
      if (inherits(object, "PNAAssay5")) {
        object <- object %>% JoinLayers()
      }
      numi <- LayerData(object, layer = "counts") %>% Matrix::colSums()
      numi <- tibble(
        component = names(numi),
        umi_count = numi %>% unname()
      )

      if (lazy) {
        copy_to(proximity_scores$src$con, numi, "numi")
        numi <- tbl(proximity_scores$src$con, "numi")
      }

      proximity_scores <- proximity_scores %>%
        left_join(numi, by = "component") %>%
        mutate(p1 = count_1 / umi_count, p2 = count_2 / umi_count) %>%
        select(-umi_count) %>%
        compute(name = "extended_proximity")
    }
  }

  return(proximity_scores)
}

#' @method ProximityScores PNAAssay5
#' @rdname ProximityScores
#'
#' @export
#'
ProximityScores.PNAAssay5 <- ProximityScores.PNAAssay


#' @method ProximityScores<- PNAAssay
#' @rdname ProximityScores
#'
#' @examples
#' # Set proximity scores
#' ProximityScores(seur_obj[["PNA"]]) <-
#'   ProximityScores(seur_obj[["PNA"]]) %>%
#'   mutate(ratio = join_count / join_count_expected_mean)
#'
#' @export
#'
"ProximityScores<-.PNAAssay" <- function(
  object,
  ...,
  value
) {
  # Validate value
  proximity <- .validate_proximity(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "proximity") <- proximity
  return(object)
}

#' @method ProximityScores<- PNAAssay5
#' @rdname ProximityScores
#'
#' @export
#'
"ProximityScores<-.PNAAssay5" <- function(
  object,
  ...,
  value
) {
  # Validate value
  proximity <- .validate_proximity(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "proximity") <- proximity
  return(object)
}


#' @param lazy A logical indicating whether to lazy load the edgelist(s)
#' from the PXL files
#'
#' @method Edgelists PNAAssay
#' @rdname Edgelists
#'
#' @examples
#' library(pixelatorR)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' el <- Edgelists(seur_obj[["PNA"]], lazy = FALSE)
#' el
#'
#' @export
#'
Edgelists.PNAAssay <- function(
  object,
  lazy = TRUE,
  ...
) {
  fs_map <- FSMap(object)
  edgelists <- .lazy_load_table(fs_map, "edgelist")

  if (!lazy) {
    con <- edgelists$src$con
    edgelists <- edgelists %>% collect()
    DBI::dbDisconnect(con)
  }

  return(edgelists)
}

#' @method Edgelists PNAAssay5
#' @rdname Edgelists
#'
#' @export
#'
Edgelists.PNAAssay5 <- Edgelists.PNAAssay

#' Get and set fs_map in a PNAAssay or PNAAssay5 object
#'
#' @param object A \code{PNAAssay} or \code{PNAAssay5} object
#' @param value A new \code{tbl_df} to replace the current fs_map with
#' @param ... Additional arguments (not used)
#'
#' @method FSMap PNAAssay
#' @rdname FSMap
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' FSMap(seur_obj[["PNA"]])
#'
#' @export
#'
FSMap.PNAAssay <- function(
  object,
  ...
) {
  slot(object, name = "fs_map")
}

#' @method FSMap PNAAssay5
#' @rdname FSMap
#'
#' @export
#'
FSMap.PNAAssay5 <- FSMap.PNAAssay


#' @method FSMap<- PNAAssay
#' @rdname FSMap
#'
#' @examples
#' # If the PXL has been moved, we can update the fs_map
#' # Here we copy the test PXL file to a temporary location
#' # to illustrate how to update the fs_map
#' temp_file <- fs::file_temp(ext = "pxl")
#' fs::file_copy(path = pxl_file, temp_file)
#'
#' # Change file path in fs_map
#' FSMap(seur_obj[["PNA"]]) <- FSMap(seur_obj[["PNA"]]) %>%
#'   mutate(pxl_file = temp_file)
#'
#' # Now the fs_map has been updated with the correct path
#' FSMap(seur_obj[["PNA"]])
#'
#' @export
#'
"FSMap<-.PNAAssay" <- function(
  object,
  ...,
  value
) {
  # Validate value
  .validate_fs_map_pna(value, data_type = "PNA")
  slot(object, name = "fs_map") <- value
  return(object)
}

#' @method FSMap<- PNAAssay5
#' @rdname FSMap
#'
#' @export
#'
"FSMap<-.PNAAssay5" <- function(
  object,
  ...,
  value
) {
  # Validate value
  .validate_fs_map_pna(value, data_type = "PNA")
  slot(object, name = "fs_map") <- value
  return(object)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' PNAAssay Methods
#'
#' Methods for \code{\link{PNAAssay}} objects for generics defined in other
#' packages
#'
#' @param x A \code{\link{PNAAssay}} object
#' @param features Feature names
#' @param cells Cell names
#' @param y A \code{\link{PNAAssay}} object or a list of \code{\link{PNAAssay}} objects
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization);
#' this is recommended if the same normalization approach was applied to all objects
#' @param ... Arguments passed to other methods
#'
#' @name PNAAssay-methods
#' @rdname PNAAssay-methods
#'
#' @concept assay
#'
NULL

#' PNAAssay5 Methods
#'
#' Methods for \code{\link{PNAAssay5}} objects for generics defined in other
#' packages
#'
#' @param x A \code{\link{PNAAssay5}} object
#' @param features Feature names
#' @param cells Cell names
#' @param y A \code{\link{PNAAssay5}} object or a list of \code{\link{PNAAssay5}} objects
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization);
#' this is recommended if the same normalization approach was applied to all objects
#' @param ... Arguments passed to other methods
#'
#' @name PNAAssay5-methods
#' @rdname PNAAssay5-methods
#'
#' @concept assay
#'
NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Base methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @describeIn PNAAssay-methods Show method for \code{PNAAssay} objects
#' @param object A \code{PNAAssay} or a \code{PNAAssay5} object
#'
#' @examples
#' library(pixelatorR)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' seur_obj[["PNA"]]
#'
#' @export
#'
setMethod(
  f = "show",
  signature = "PNAAssay",
  definition = function(object) {
    cellgraphs <- slot(object, "cellgraphs")
    loaded_graphs <- !sapply(cellgraphs, is.null)
    msg <-
      capture.output(show(as(
        object,
        Class = ifelse(
          inherits(object, "PNAAssay"),
          "Assay",
          "Assay5"
        )
      )))
    msg[1] <- gsub(
      pattern = "^Assay",
      x = msg[1],
      replacement = "PNAAssay"
    )
    msg <- c(msg, paste0(
      "Loaded CellGraph objects:\n ",
      sum(loaded_graphs),
      "\n"
    ))
    cat(paste(msg, collapse = "\n"))
  }
)


#' @describeIn PNAAssay5-methods Show method for \code{PNAAssay5} objects
#' @concept assay
#' @method show CellGraphAssay
#' @docType methods
#'
#' @export
#'
setMethod(
  f = "show",
  signature = "PNAAssay5",
  definition = function(object) as(getMethod("show", "PNAAssay"), "function")(object)
)


#' @describeIn PNAAssay-methods Subset a \code{PNAAssay} or a
#' \code{PNAAssay5} object
#' @concept assay
#' @method subset PNAAssay
#' @docType methods
#'
#' @return A \code{PNAAssay} object
#'
#' @examples
#' library(pixelatorR)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' pna_assay <- seur_obj[["PNA"]]
#'
#' pna_assay <- subset(pna_assay, cells = colnames(pna_assay)[1:2])
#' pna_assay
#'
#' @export
#'
subset.PNAAssay <- function(
  x,
  features = NULL,
  cells = NULL,
  ...
) {
  assert_x_in_y(cells, colnames(x), allow_null = TRUE)
  assert_x_in_y(features, rownames(x), allow_null = TRUE)

  # Get cellgraphs
  cellgraphs <- x@cellgraphs

  # subset elements in the standard assay
  assay <- as(object = x, Class = ifelse(is(x, "PNAAssay"),
    "Assay",
    "Assay5"
  ))
  if (inherits(assay, "Assay5")) {
    assay <- JoinLayers(assay)
  }
  assay_subset <- subset(x = assay, features = features, cells = cells)

  # Filter cellgraphs
  cellgraphs_filtered <- cellgraphs[colnames(assay_subset)]

  # Filter proximity scores
  proximity <- slot(x, name = "proximity")
  if (length(proximity) > 0) {
    if (!is.null(cells)) {
      proximity <- proximity %>% filter(component %in% cells)
    }
    if (!is.null(features)) {
      proximity <- proximity %>% filter(marker_1 %in% features & marker_2 %in% features)
    }
  }

  # Filter fs_map
  fs_map <- slot(x, name = "fs_map")
  if (length(fs_map) > 0) {
    fs_map$id_map <- lapply(fs_map$id_map, function(x) {
      x %>% filter(current_id %in% cells)
    })
    # Drop sample if id_map is empty
    rows_keep <- sapply(fs_map$id_map, nrow) > 0
    fs_map <- fs_map[rows_keep, ]
    fs_map <- na.omit(fs_map)
    # Update sample integer
    fs_map$sample <- seq_len(nrow(fs_map))
  }

  # convert standard assay to PNAAssay or PNAAssayy5
  as_assay_func <- ifelse(is(x, "PNAAssay"),
    as.PNAAssay,
    as.PNAAssay5
  )
  pna_assay <- as_assay_func(
    x = assay_subset,
    cellgraphs = cellgraphs_filtered,
    proximity = proximity,
    fs_map = fs_map
  )
  return(pna_assay)
}


#' @describeIn PNAAssay5-methods Subset a \code{PNAAssay5} object
#' @concept assay
#' @method subset PNAAssay5
#' @docType methods
#'
#' @export
#'
subset.PNAAssay5 <- subset.PNAAssay


#' @param add.cell.ids A character vector with sample names
#' @param collapse If TRUE, merge layers of the same name together
#'
#' @describeIn PNAAssay-methods Merge two or more \code{PNAAssay} or
#' \code{PNAAssay5} objects together
#' @concept assay
#' @method merge PNAAssay
#' @docType methods
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- minimal_pna_pxl_file()
#' seur_obj <- ReadPNA_Seurat(pxl_file)
#' pna_assay <- seur_obj[["PNA"]]
#'
#' # Merge two data sets
#' pna_assay_merged <-
#'   merge(pna_assay, pna_assay, add.cell.ids = c("A", "B"))
#' pna_assay_merged
#'
#' # Merge multiple data sets
#' pna_assay_list <- list(pna_assay, pna_assay, pna_assay)
#' pna_assay_merged <-
#'   merge(pna_assay_list[[1]], pna_assay_list[-1], add.cell.ids = c("A", "B", "C"))
#' pna_assay_merged
#'
#' @export
#'
merge.PNAAssay <- function(
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  collapse = TRUE,
  ...
) {
  # Validate input parameters
  assert_class(y, c("list", "PNAAssay", "PNAAssay5"))
  if (is.list(y)) {
    for (i in seq_along(y)) {
      if (!inherits(y[[i]], what = c("PNAAssay", "PNAAssay5"))) {
        cli::cli_abort(
          c(
            "i" = "All elements of {.var y} must be {.cls CellGraphAssay} or {.cls CellGraphAssay5} objects",
            "x" = "Element {i} of {.var y} is a {.cls {class(y[[i]])}}"
          )
        )
      }
    }
  }

  objects <- c(x, y)
  cell_names <- unlist(lapply(objects, colnames))

  # Define add.cell.ids
  if (!is.null(add.cell.ids)) {
    if (!(length(add.cell.ids) == length(objects))) {
      cli::cli_abort(
        c(
          "i" = "Length of 'add.cell.ids' must match the number of objects to merge",
          "x" = "Length of 'add.cell.ids': {.val {length(add.cell.ids)}}",
          "x" = "Number of objects to merge: {.val {length(objects)}}"
        )
      )
    }
    objects <- lapply(seq_along(objects), function(i) {
      objects[[i]] %>%
        RenameCells(new.names = paste0(add.cell.ids[i], "_", colnames(objects[[i]])) %>%
          set_names(nm = colnames(objects[[i]])))
    })
  }

  # Check duplicate cell names
  unique_names <- table(cell_names)
  names_are_duplicated <- any(unique_names > 1)
  if (names_are_duplicated && is.null(add.cell.ids)) {
    cli::cli_abort(
      c(
        "i" = "Found non-unique cell IDs across samples. {.var add.cell.ids} must be specified. ",
        " " = "Alternatively, make sure that each separate object has unique cell IDs.",
        "x" = "Duplicated IDs: {.val {names(unique_names)[unique_names > 1]}}"
      )
    )
  }

  # Fetch standard assays
  standardassays <- lapply(objects, function(pna_assay) {
    assay <- as(object = pna_assay, Class = ifelse(is(pna_assay, "PNAAssay"), "Assay", "Assay5"))
    if (length(Key(assay)) == 0) Key(assay) <- "px_"
    return(assay)
  })

  # Merge Seurat assays
  new_assay <- merge(
    x = standardassays[[1]],
    y = standardassays[-1],
    merge.data = merge.data,
    ...
  )

  # Join layers
  if (collapse && is(new_assay, "PNAAssay5")) {
    new_assay <- JoinLayers(new_assay)
  }
  if (collapse && is(new_assay, "Assay5")) {
    new_assay <- JoinLayers(new_assay)
  }

  # Merge cellgraphs list
  cellgraphs_new <- Reduce(c, lapply(objects, function(pna_assay) {
    return(pna_assay@cellgraphs)
  })) %>% set_names(nm = colnames(new_assay))

  # Merge spatial metrics tables
  proximity <- lapply(objects, function(object) {
    slot(object, name = "proximity")
  }) %>% bind_rows()

  # Merge fs_map
  fs_map <- tibble()
  for (pna_assay in objects) {
    if (length(slot(pna_assay, name = "fs_map")) == 0) {
      fs_map <- NULL
      break
    }
  }
  if (!is.null(fs_map)) {
    fs_map <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
      slot(objects[[i]], name = "fs_map")
    })) %>% mutate(sample = seq_len(n()))
  }

  # Convert to PNAAssay and add cellgraphs
  as_assay_func <-
    ifelse(is(pna_assay, "PNAAssay"),
      as.PNAAssay,
      as.PNAAssay5
    )
  merged_pna_assay <- as_assay_func(
    x = new_assay,
    cellgraphs = cellgraphs_new,
    proximity = proximity,
    fs_map = fs_map
  )

  # Add meta.features for PNAAssay
  if (is(merged_pna_assay, "PNAAssay")) {
    can_merge_meta_features <- TRUE
    ob_features <- lapply(objects, function(x) rownames(x@meta.features))
    for (i in 2:length(ob_features)) {
      if (!all(ob_features[[i]] == ob_features[[1]])) {
        cli_alert_danger(
          "Meta features are not the same across objects. ",
          "Cannot merge meta.features slots."
        )
        can_merge_meta_features <- FALSE
        break
      }
    }

    if (can_merge_meta_features) {
      meta_features <- objects[[1]]@meta.features %>% select(any_of(c("control", "marker", "nuclear")))
    }

    # Place meta.features in CellGraphAssay
    merged_pna_assay@meta.features <- meta_features
  }

  return(merged_pna_assay)
}

#' @describeIn PNAAssay5-methods Merge two or more \code{PNAAssay5} objects together
#' @concept assay
#' @method merge PNAAssay5
#' @docType methods
#'
#' @export
#'
merge.PNAAssay5 <- merge.PNAAssay


#' Validate proximity \code{tbl_df}
#'
#' @param proximity A \code{tbl_df} with proximity scores
#' @param cell_ids A character vector with cell IDs
#' @param markers A character vector with marker names
#' @param call Environment to use for the error call
#' @param verbose Print messages
#'
#' @noRd
.validate_proximity <- function(
  proximity,
  cell_ids,
  markers,
  call = caller_env(),
  verbose = FALSE
) {
  # Set proximity to empty tibble if NULL
  proximity <- proximity %||% tibble()
  assert_class(proximity, "tbl_df", call = call)

  if (length(proximity) > 0) {
    # Check column names
    name_check <- all(c("marker_1", "marker_2", "component") %in% names(proximity))
    if (!name_check) {
      cli::cli_abort(
        c(
          "i" = "Columns {.str marker_1}, {.str marker_2} and {.str component} must",
          " " = "be present in the 'proximity' score table"
        ),
        call = call
      )
    }
    # Check component names
    cells_in_proximity <- cell_ids %in% (proximity$component %>% unique())
    if (!all(cells_in_proximity)) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(glue(
          "Cells {paste(which(!cells_in_proximity), collapse=', ')} ",
          "in count matrix are missing from {.str proximity} table"
        ))
      }
      return(proximity)
    }
    # Check marker names
    all_markers <- c(proximity$marker_1 %>% unique(), proximity$marker_2 %>% unique()) %>% unique()
    if (!all(markers %in% all_markers)) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(glue(
          "Markers {paste(which(!markers %in% all_markers), collapse=', ')} ",
          "in count matrix are missing from 'proximity' table"
        ))
      }
    }
  }
  return(proximity)
}


#' Validate an fs_map tibble
#'
#' TODO: update .validate_fs_map to this function in pixelatorR when merging
#'
#' @noRd
.validate_fs_map_pna <- function(
  fs_map,
  allow_empty_df = FALSE,
  data_type = c("MPX", "PNA"),
  call = caller_env()
) {
  if (inherits(fs_map, "data.frame") && length(fs_map) == 0 && allow_empty_df) {
    return(invisible(NULL))
  }
  data_type <- match.arg(data_type, choices = c("MPX", "PNA"))
  assert_non_empty_object(fs_map, "tbl_df", call = call)
  if (!all(c("id_map", "sample", "pxl_file") == colnames(fs_map))) {
    cli::cli_abort(
      c(
        "i" = "{.var fs_map} must have columns {.str id_map}, {.str sample}, and {.str pxl_file}",
        "x" = "Found column(s) {.val {colnames(fs_map)}} in {.var fs_map}"
      ),
      call = call
    )
  }

  # Validate columns
  assert_col_class("id_map", fs_map, classes = "list", call = call)
  assert_col_class("sample", fs_map, classes = "integer", call = call)
  assert_col_class("pxl_file", fs_map, classes = "character", call = call)
  id_map_check1 <- sapply(fs_map$id_map, function(x) inherits(x, what = "tbl_df"))
  if (!all(id_map_check1)) {
    cli::cli_abort(
      c("x" = "All elements of {.var id_map} must be {.cls tbl_df} objects"),
      call = call
    )
  }
  id_map_check2 <- sapply(fs_map$id_map, function(x) all(colnames(x) == c("current_id", "original_id")))
  if (!all(id_map_check2)) {
    cli::cli_abort(
      c(
        "x" = "All elements of {.var id_map} must be {.cls tbl_df} objects with",
        " " = "columns {.str current_id} and {.str original_id}"
      ),
      call = call
    )
  }
  id_map_check3 <- sapply(fs_map$id_map, function(x) {
    all(sapply(x, class) == c("character", "character"))
  })
  if (!all(id_map_check3)) {
    cli::cli_abort(
      c("x" = "Invalid classes found in {.var id_map} columns 'current_id' and 'original_id'"),
      call = call
    )
  }

  # Validate pxl files
  for (i in seq_len(nrow(fs_map))) {
    f <- fs_map$pxl_file[i]
    if (!fs::file_exists(f)) {
      cli::cli_abort(
        c(
          "x" = "The pxl file {.file {f}} linked to sample {.val {i}} does not exist. ",
          "i" = "Make sure that the path is correct and that the file has not been moved/deleted."
        ),
        call = call
      )
    }
  }

  return(invisible(NULL))
}
