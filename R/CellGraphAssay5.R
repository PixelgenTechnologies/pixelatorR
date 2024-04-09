#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod slot slot<- new as slotNames
#' @importClassesFrom Matrix dgCMatrix
NULL

# Declarations used in package check
globalVariables(
  names = c('marker_1', 'marker_2', 'id_map', 'component_new',
            'current_id'),
  package = 'pixelatorR',
  add = TRUE
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The CellGraphAssay5 class
#'
#' The CellGraphAssay5 object is an extended \code{\link[SeuratObject]{Assay5}}
#' for the storage and analysis of mpx single-cell data.
#'
#' @slot cellgraphs A named list of \code{\link{CellGraph}} objects
#' @slot polarization A \code{tbl_df} with polarization scores
#' @slot colocalization A \code{tbl_df} with colocalization scores
#' @slot fs_map A \code{tbl_df} with information pxl file paths,
#' sample IDs and component IDs
#'
#' @name CellGraphAssay5-class
#' @rdname CellGraphAssay5-class
#' @importClassesFrom SeuratObject Assay5
#' @exportClass CellGraphAssay5
#' @concept assay
CellGraphAssay5 <- setClass(
  Class = "CellGraphAssay5",
  contains = "Assay5",
  slots = list(
    cellgraphs = "list",
    polarization = "data.frame",
    colocalization = "data.frame",
    fs_map = "data.frame"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create a CellGraphAssay5 object
#'
#' Create a \code{\link{CellGraphAssay5}} object from a count matrix. The expected
#' format of the input matrix is features x cells.
#'
#' @param counts Unnormalized data (raw counts)
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects
#' @param polarization A \code{tbl_df} with polarization scores
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param fs_map A \code{tbl_df} with information pxl file paths,
#' sample IDs and component IDs
#' @param ... Additional arguments passed to \code{\link{CreateAssay5Object}}
#' @inheritParams ReadMPX_arrow_edgelist
#'
#' @import rlang
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom Matrix rowSums colSums
#' @concept assay
#'
#' @return A \code{CellGraphAssay5} object
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' counts <- ReadMPX_counts(pxl_file)
#' edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#' components <- colnames(counts)
#' edgelist_split <-
#'   edgelist %>%
#'   select(upia, upib, component) %>%
#'   distinct() %>%
#'   group_by(component) %>%
#'   group_split() %>%
#'   setNames(nm = components)
#'
#' # Convert data into a list of CellGraph objects
#' bipartite_graphs <- lapply(edgelist_split, function(x) {
#'   x <- x %>% as_tbl_graph(directed = FALSE)
#'   x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
#'   attr(x, "type") <- "bipartite"
#'   CreateCellGraphObject(cellgraph = x)
#' })
#'
#' # Create CellGraphAssay5
#' cg_assay5 <- CreateCellGraphAssay5(counts = counts, cellgraphs = bipartite_graphs)
#' cg_assay5
#'
#' @export
#'
CreateCellGraphAssay5 <- function (
  counts,
  cellgraphs,
  polarization = NULL,
  colocalization = NULL,
  fs_map = NULL,
  verbose = FALSE,
  ...
) {

  # Check input parameters
  stopifnot(
    "'counts' must be a matrix-like object" =
      inherits(counts, what =  c("matrix", "dgCMatrix")),
    "'cellgraphs' must be a 'list'" =
      inherits(cellgraphs, what = "list")
  )
  stopifnot(
    "'cellgraphs' names must be match the colnames of the count matrix" =
      all(names(cellgraphs) == colnames(counts))
  )
  cellgraphs <- cellgraphs[colnames(counts)]

  # Validate polarization and colocalization
  polarization <-
    .validate_polarization(
      polarization,
      cell_ids = colnames(counts),
      markers = rownames(counts),
      verbose = verbose
    )
  colocalization <-
    .validate_colocalization(
      colocalization,
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
  cg_assay5 <- as.CellGraphAssay5(
    x = seurat_assay,
    cellgraphs = cellgraphs,
    polarization = polarization,
    colocalization = colocalization,
    fs_map = fs_map
  )
  return(cg_assay5)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get/set methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @rdname CellGraphs
#' @method CellGraphs CellGraphAssay5
#' @export
#' @concept assay
#' @concept cellgraphs
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' counts <- ReadMPX_counts(pxl_file)
#' edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#' components <- colnames(counts)
#' edgelist_split <-
#'   edgelist %>%
#'   select(upia, upib, component) %>%
#'   distinct() %>%
#'   group_by(component) %>%
#'   group_split() %>%
#'   setNames(nm = components)
#'
#' # Convert data into a list of CellGraph objects
#' bipartite_graphs <- lapply(edgelist_split, function(x) {
#'   x <- x %>% as_tbl_graph(directed = FALSE)
#'   x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
#'   attr(x, "type") <- "bipartite"
#'   CreateCellGraphObject(cellgraph = x)
#' })
#'
#' # CellGraphs getter CellGraphAssay5
#' # ---------------------------------
#'
#' # Create CellGraphAssay5
#' cg_assay5 <- CreateCellGraphAssay5(counts = counts, cellgraphs = bipartite_graphs)
#' cg_assay5
#'
#' # Get cellgraphs from a CellGraphAssay5 object
#' CellGraphs(cg_assay5)
#'
CellGraphs.CellGraphAssay5 <- function (
  object,
  ...
) {
  return(slot(object, name = "cellgraphs"))
}


#' @export
#' @method CellGraphs<- CellGraphAssay5
#' @rdname CellGraphs
#' @concept assay
#' @concept cellgraphs
#'
#' @examples
#'
#' # CellGraphs setter for CellGraphAssay5
#' # --------------------------------------
#'
#' # Set cellgraphs in a CellGraphAssay5 object
#' CellGraphs(cg_assay5) <- cg_assay5@cellgraphs
#'
"CellGraphs<-.CellGraphAssay5" <- function (
  object,
  ...,
  value
) {

  # Clean cellgraphs slot if balue = NULL
  if (is.null(x = value)) {
    slot(object = object, name = "cellgraphs") <- rep(list(NULL), ncol(object)) %>% setNames(nm = colnames(object))
    return(object)
  }

  # Validate list
  if (inherits(x = value, what = "list")) {
    if (!all(names(value) == names(CellGraphs(object))))
      abort("new list names must match the current cellgraphs list names")
    for (i in seq_along(value)) {
      if (!inherits(x = value[[i]], what = c("CellGraph", "NULL"))) {
        abort(glue("Element {i} is not a CellGraph object or NULL"))
      }
    }
    slot(object = object, name = "cellgraphs") <- value
  } else {
    abort(glue("Invalid class '{class(value)}'"))
  }
  return(object)
}


#' @importFrom SeuratObject RenameCells LayerData
#'
#' @rdname RenameCells
#' @method RenameCells CellGraphAssay5
#' @concept assay
#'
#' @export
RenameCells.CellGraphAssay5 <- function (
  object,
  new.names,
  ...
) {

  if (!inherits(new.names, what = "character") ||
      !(length(new.names) == ncol(object)) ||
      !sum(duplicated(new.names)) == 0) {
    abort(glue("'new.names' must be a character vector where length(new.names) == ncol(object),",
               " and the names must be unique."))
  }
  if (is.null(names(new.names))) {
    new.names <- set_names(new.names, nm = colnames(object))
  }

  orig.names <- colnames(object)

  assay5 <- as(object, Class = "Assay5")
  assay5 <- RenameCells(assay5, new.names = new.names, ...)

  ## Recreate the CellGraphAssay5 object
  cg_assay5_renamed <- as.CellGraphAssay5(assay5)

  # Rename cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  cellgraphs <- set_names(cellgraphs, nm = new.names)
  slot(cg_assay5_renamed, name = "cellgraphs") <- cellgraphs

  # Handle polarization slot
  name_conversion <- tibble(new = new.names, component = orig.names)
  polarization <- slot(object, name = "polarization")
  if (length(polarization) > 0) {
    polarization <- polarization %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(cg_assay5_renamed, name = "polarization") <- polarization

  # Handle colocalization slot
  colocalization <- slot(object, name = "colocalization")
  if (length(colocalization) > 0) {
    colocalization <- colocalization %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(cg_assay5_renamed, name = "colocalization") <- colocalization

  # Handle fs_map
  fs_map <- slot(object, name = "fs_map")
  if (nrow(fs_map) > 0) {
    fs_map$id_map <- fs_map %>% pull(id_map) %>%
      lapply(function(x) {
        x$current_id <- new.names[x$current_id]
        return(x)
      })
  }
  slot(cg_assay5_renamed, name = "fs_map") <- fs_map

  return(cg_assay5_renamed)
}

#' @param cellgraphs A list of \code{\link{CellGraph}} objects
#' @param polarization A \code{tbl_df} with polarization scores
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param fs_map A \code{tbl_df} with information pxl file paths,
#' sample IDs and component IDs
#'
#' @import rlang
#'
#' @rdname as.CellGraphAssay5
#' @method as.CellGraphAssay5 Assay5
#' @concept assay
#'
#' @return A \code{CellGraphAssay5} object
#'
#' @examples
#'
#' library(pixelatorR)
#' library(SeuratObject)
#' library(dplyr)
#' library(tidygraph)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' counts <- ReadMPX_counts(pxl_file)
#' edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#' components <- colnames(counts)
#' edgelist_split <-
#'   edgelist %>%
#'   select(upia, upib, component) %>%
#'   distinct() %>%
#'   group_by(component) %>%
#'   group_split() %>%
#'   setNames(nm = components)
#'
#' # Convert data into a list of CellGraph objects
#' bipartite_graphs <- lapply(edgelist_split, function(x) {
#'   x <- x %>% as_tbl_graph(directed = FALSE)
#'   x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
#'   attr(x, "type") <- "bipartite"
#'   CreateCellGraphObject(cellgraph = x)
#' })
#'
#' # Create Assay5
#' assay5 <- CreateAssay5Object(counts = counts)
#'
#' # Convert Assay5 to CellGraphAssay5
#' cg_assay <- as.CellGraphAssay5(assay5, cellgraphs = bipartite_graphs)
#' cg_assay
#'
#' @export
#'
as.CellGraphAssay5.Assay5 <- function (
  x,
  cellgraphs = NULL,
  polarization = NULL,
  colocalization = NULL,
  fs_map = NULL,
  ...
) {

  # Check cellgraphs
  if (!is.null(cellgraphs)) {
    stopifnot(
      "'cellgraphs' must be a non-empty list with the same number of elements as the number of columns in the Assay5" =
        is.list(cellgraphs) &&
        (length(cellgraphs) == ncol(x)))
    stopifnot(
      "'cellgraphs' names must match colnames of the Assay5" =
        all(names(cellgraphs) == colnames(x))
    )
    for (i in seq_along(cellgraphs)) {
      if (!inherits(x = cellgraphs[[i]], what = c("CellGraph", "NULL"))) {
        abort(glue("Element {i} is not a CellGraph object or NULL"))
      }
    }
  } else {
    cellgraphs <- rep(list(NULL), ncol(x)) %>% setNames(nm = colnames(x))
  }

  # Check fs_map
  if (!is.null(fs_map)) {
    .validate_fs_map(fs_map)
  }

  # Validate polarization and colocalization
  polarization <- .validate_polarization(polarization, cell_ids = colnames(x), markers = rownames(x))
  colocalization <- .validate_colocalization(colocalization, cell_ids = colnames(x), markers = rownames(x))

  new_assay <- as(object = x, Class = "CellGraphAssay5")

  # Add slots
  slot(new_assay, name = "cellgraphs") <- cellgraphs
  slot(new_assay, name = "polarization") <- polarization
  slot(new_assay, name = "colocalization") <- colocalization

  # Add fs_map
  if (!is.null(fs_map)) {
    slot(new_assay, name = "fs_map") <- fs_map
  }

  return(new_assay)
}

setAs(
  from = "Assay5",
  to = "CellGraphAssay5",
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
        "Class" = "CellGraphAssay5"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)


#' @method PolarizationScores CellGraphAssay5
#'
#' @rdname PolarizationScores
#'
#' @export
#'
PolarizationScores.CellGraphAssay5 <- function (
  object,
  ...
) {
  slot(object, name = "polarization")
}


#' @method PolarizationScores<- CellGraphAssay5
#'
#' @rdname PolarizationScores
#'
#' @export
#'
"PolarizationScores<-.CellGraphAssay5" <- function (
  object,
  ...,
  value
) {
  # Validate value
  polarization <- .validate_polarization(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "polarization") <- polarization
  return(object)
}


#' @method ColocalizationScores CellGraphAssay5
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
ColocalizationScores.CellGraphAssay5 <- function (
  object,
  ...
) {
  slot(object, name = "colocalization")
}


#' @method ColocalizationScores<- CellGraphAssay5
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
"ColocalizationScores<-.CellGraphAssay5" <- function (
  object,
  ...,
  value
) {
  # Validate value
  colocalization <- .validate_colocalization(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "colocalization") <- colocalization
  return(object)
}


#' @method FSMap CellGraphAssay5
#'
#' @rdname FSMap
#'
#' @export
#'
FSMap.CellGraphAssay5 <- function (
  object,
  ...
) {
  slot(object, name = "fs_map")
}


#' @method FSMap<- CellGraphAssay5
#'
#' @rdname FSMap
#'
#' @export
#'
"FSMap<-.CellGraphAssay5" <- function (
  object,
  ...,
  value
) {
  # Validate value
  .validate_fs_map(value)
  slot(object, name = "fs_map") <- value
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' CellGraphAssay5 Methods
#'
#' Methods for \code{\link{CellGraphAssay5}} objects for generics defined in other
#' packages
#'
#' @param x A \code{\link{CellGraphAssay5}} object
#' @param features Feature names
#' @param cells Cell names
#' @param y A \code{\link{CellGraphAssay5}} object or a list of \code{\link{CellGraphAssay5}} objects
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization);
#' this is recommended if the same normalization approach was applied to all objects
#' @param ... Arguments passed to other methods
#'
#' @name CellGraphAssay5-methods
#' @rdname CellGraphAssay5-methods
#'
#' @concept assay
#'
NULL


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Base methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Show method for \code{CellGraphAssay5} object
#'
#' @param object A \code{CellGraphAssay5} object
#'
#' @importFrom methods show
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' counts <- ReadMPX_counts(pxl_file)
#' edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#' components <- colnames(counts)
#' edgelist_split <-
#'   edgelist %>%
#'   select(upia, upib, component) %>%
#'   distinct() %>%
#'   group_by(component) %>%
#'   group_split() %>%
#'   setNames(nm = components)
#'
#' # Convert data into a list of CellGraph objects
#' bipartite_graphs <- lapply(edgelist_split, function(x) {
#'   x <- x %>% as_tbl_graph(directed = FALSE)
#'   x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
#'   attr(x, "type") <- "bipartite"
#'   CreateCellGraphObject(cellgraph = x)
#' })
#'
#' # Create CellGraphAssay5
#' cg_assay5 <- CreateCellGraphAssay5(counts = counts, cellgraphs = bipartite_graphs)
#' cg_assay5
#'
#' # Show method
#' cg_assay5
#'
setMethod (
  f = "show",
  signature = "CellGraphAssay5",
  definition = function(object) {
    cellgraphs <- slot(object, "cellgraphs")
    loaded_graphs <- !sapply(cellgraphs, is.null)
    show(as(object, Class = "Assay5"))
    cat(
      "Loaded CellGraph objects:",
      sum(loaded_graphs),
      "\n"
    )
  }
)



#' @describeIn CellGraphAssay5-methods Subset a \code{CellGraphAssay5} object
#' @importClassesFrom SeuratObject Assay
#' @concept assay
#' @method subset CellGraphAssay5
#'
#' @importFrom arrow write_dataset
#'
#' @return A \code{CellGraphAssay5} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' counts <- ReadMPX_counts(pxl_file)
#' edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#' components <- colnames(counts)
#' edgelist_split <-
#'   edgelist %>%
#'   select(upia, upib, component) %>%
#'   distinct() %>%
#'   group_by(component) %>%
#'   group_split() %>%
#'   setNames(nm = components)
#'
#' # Convert data into a list of CellGraph objects
#' bipartite_graphs <- lapply(edgelist_split, function(x) {
#'   x <- x %>% as_tbl_graph(directed = FALSE)
#'   x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
#'   attr(x, "type") <- "bipartite"
#'   CreateCellGraphObject(cellgraph = x)
#' })
#'
#' # Create CellGraphAssay5
#' cg_assay5 <- CreateCellGraphAssay5(counts = counts,
#'                                    cellgraphs = bipartite_graphs,
#'                                    fs_map = tibble(
#'                                      id_map = list(tibble(
#'                                        current_id = colnames(counts),
#'                                        original_id = colnames(counts)
#'                                      )), sample = 1L, pxl_file = pxl_file
#'                                    ))
#' cg_assay5 <- LoadCellGraphs(cg_assay5)
#'
#' # Subset CellGraphAssay5
#' # ---------------------------------
#' cg_assay_subset <- subset(cg_assay5, cells = colnames(cg_assay5)[1:3])
#'
#' @export
#'
subset.CellGraphAssay5 <- function (
  x,
  features = NULL,
  cells = NULL,
  ...
) {

  stopifnot(
    "All 'cells' must be present in x" =
      all(cells %in% colnames(x))
  )

  # Get cellgraphs
  cellgraphs <- x@cellgraphs

  # subset elements in the standard assay
  assay5 <- as(object = x, Class = "Assay5")
  assay5_subset <- subset(x = assay5,
                   features = features,
                   cells = cells)

  # Filter cellgraphs
  cellgraphs_filtered <- cellgraphs[colnames(assay5_subset)]

  # Filter polarization and colocaliation scores
  polarization <- slot(x, name = "polarization")
  if (length(polarization) > 0) {
    if (!is.null(cells)) {
      polarization <- polarization %>% filter(component %in% cells)
    }
    if (!is.null(features)) {
      polarization <- polarization %>% filter(marker %in% features)
    }
  }
  colocalization <- slot(x, name = "colocalization")
  if (length(colocalization) > 0) {
    if (!is.null(cells)) {
      colocalization <- colocalization %>% filter(component %in% cells)
    }
    if (!is.null(features)) {
      colocalization <- colocalization %>% filter((marker_1 %in% features) & (marker_2 %in% features))
    }
  }

  # Filter fs_map
  fs_map <- slot(x, name = "fs_map")
  if (length(fs_map) > 0) {
    fs_map$id_map <- lapply(fs_map$id_map, function(x) {
      x %>% filter(current_id %in% cells)
    })
    fs_map <- na.omit(fs_map)
  }

  # convert standard assay to CellGraphAssay5
  cg_assay5 <- as.CellGraphAssay5(
    x = assay5_subset,
    cellgraphs = cellgraphs_filtered,
    polarization = polarization,
    colocalization = colocalization,
    fs_map = fs_map
  )
  return(cg_assay5)
}


#' @param add.cell.ids A character vector with sample names
#'
#' @importFrom SeuratObject Key Key<- RenameCells
#' @describeIn CellGraphAssay5-methods Merge two or more \code{CellGraphAssay5} objects together
#' @concept assay
#' @method merge CellGraphAssay5
#'
#' @importFrom SeuratObject Key Key<- RenameCells
#'
#' @examples
#' # Merge CellGraphAssay5s
#' # ---------------------------------
#'
#' # Merge 3 CellGraphAssay5s
#' cg_assay5_merged <- merge(cg_assay5,
#'                           y = list(cg_assay5, cg_assay5),
#'                           add.cell.ids = c("A", "B", "C"))
#'
#' @export
#'
merge.CellGraphAssay5 <- function (
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  ...
) {

  # Validate input parameters
  if (!inherits(y, what = c("list", "CellGraphAssay5")))
    abort("'y' must be a 'CellGraphAssay5' object or a list of 'CellGraphAssay5' objects")
  if (is.list(y)) {
    for (i in seq_along(y)) {
      if (!inherits(y[[i]], what = "CellGraphAssay5")) {
        abort(glue("Element {i} in 'y' is not a 'CellGraphAssay5'"))
      }
    }
  }

  objects <- c(x, y)
  cell.names <- unlist(lapply(objects, colnames))
  name_conversion <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
    tibble(component = colnames(objects[[i]]), sample = i)
  })) %>% mutate(component_new = cell.names) %>%
    group_by(sample) %>%
    group_split()

  # Define add.cell.ids
  if (!is.null(add.cell.ids)) {
    stopifnot(
      "Length of 'add.cell.ids' must match the number of objects to merge" =
        length(add.cell.ids) == length(objects)
    )
    objects <- lapply(seq_along(objects), function(i) {
      objects[[i]] %>% RenameCells(new.names = paste0(add.cell.ids[i], "_", colnames(objects[[i]])) %>%
                                     set_names(nm = colnames(objects[[i]])))
    })
  }

  # Check duplicate cell names
  unique_names <- table(cell.names)
  names_are_duplicated <- any(unique_names > 1)
  if (names_are_duplicated & is.null(add.cell.ids)) {
    abort(glue("Found non-unique IDs across samples. A 'add.cell.ids' must be specified. ",
               "Alternatively, make sure that each separate object has unique cell names."))
  }

  # Fetch standard assays
  standardassays <- lapply(objects, function(cg_assay5) {
    assay5 <- as(object = cg_assay5, Class = "Assay5")
    if (length(Key(assay5)) == 0) Key(assay5) <- "px_"
    return(assay5)
  })

  # Merge Seurat assays
  new_assay5 <- merge(x = standardassays[[1]],
                      y = standardassays[-1],
                      merge.data = merge.data,
                      ...)

  # Fetch cellgraphs
  cellgraphs_new <- Reduce(c, lapply(objects, function(cg_assay5) {
    return(cg_assay5@cellgraphs)
  })) %>% set_names(nm = colnames(new_assay5))


  # Merge polarization and colocalization scores
  polarization <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
    pl <- slot(objects[[i]], name = "polarization")
    if (length(pl) == 0) return(pl)
    pl <- pl %>% left_join(name_conversion[[i]], by = "component") %>%
      select(-component, -sample) %>%
      rename(component = component_new)
    return(pl)
  }))
  colocalization <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
    cl <- slot(objects[[i]], name = "colocalization")
    if (length(cl) == 0) return(cl)
    cl <- cl %>% left_join(name_conversion[[i]], by = "component") %>%
      select(-component, -sample) %>%
      rename(component = component_new)
    return(cl)
  }))

  # Merge fs_map
  fs_map <- tibble()
  for (cg_assay5 in objects) {
    if (length(slot(cg_assay5, name = "fs_map")) == 0) {
      fs_map <- NULL
      break
    }
  }
  if (!is.null(fs_map)) {
    fs_map <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
      slot(objects[[i]], name = "fs_map")
    })) %>% mutate(sample = seq_len(n()))
  }

  # Convert to CellGraphAssay5 and add cellgraphs
  merged_cg_assay5 <- as.CellGraphAssay5(
    x = new_assay5,
    cellgraphs = cellgraphs_new,
    polarization = polarization,
    colocalization = colocalization,
    fs_map = fs_map
  )

  return(merged_cg_assay5)
}
