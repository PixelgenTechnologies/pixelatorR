#' @include generics.R
#' @importClassesFrom Matrix dgCMatrix
NULL

# Declarations used in package check
globalVariables(
  names = c('marker_1', 'marker_2'),
  package = 'pixelatorR',
  add = TRUE
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The CellGraphAssay class
#'
#' The CellGraphAssay object is an extended \code{\link[SeuratObject]{Assay}}
#' for the storage and analysis of MPX single-cell data.
#'
#' @slot cellgraphs A named list of \code{\link{CellGraph}} objects
#' @slot polarization A \code{tbl_df} with polarization scores
#' @slot colocalization A \code{tbl_df} with colocalization scores
#' @slot fs_map A \code{tbl_df} with information pxl file paths,
#' sample IDs and component IDs
#'
#' @name CellGraphAssay-class
#' @rdname CellGraphAssay-class
#' @importClassesFrom SeuratObject Assay
#' @exportClass CellGraphAssay
#' @concept assay
CellGraphAssay <- setClass(
  Class = "CellGraphAssay",
  contains = "Assay",
  slots = list(
    cellgraphs = "list",
    polarization = "data.frame",
    colocalization = "data.frame",
    fs_map = "data.frame"
  )
)

#' The CellGraphAssay5 class
#'
#' The CellGraphAssay5 object is an extended \code{\link[SeuratObject]{Assay5}}
#' for the storage and analysis of MPX single-cell data.
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


#' MPXAssay class
#'
#' The MPXAssay class is defined as the union of the CellGraphAssay and CellGraphAssay5
#' class. In other words, it is a virtual class defined as a superclass of these two classes.
#'
#' @name MPXAssay-class
#' @rdname MPXAssay-class
#' @aliases MPXAssay
#'
#' @exportClass MPXAssay
#' @concept assay
setClassUnion("MPXAssay", c("CellGraphAssay", "CellGraphAssay5"))



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Create a CellGraphAssay object
#'
#' Create a \code{\link{CellGraphAssay}} object from a count matrix. The expected
#' format of the input matrix is features x cells.
#'
#' @param counts Unnormalized data (raw counts)
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects
#' @param polarization A \code{tbl_df} with polarization scores
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param fs_map A \code{tbl_df} with information pxl file paths,
#' sample IDs and component IDs
#' @param ... Additional arguments passed to \code{\link{CreateAssayObject}}
#' @inheritParams ReadMPX_arrow_edgelist
#'
#' @import rlang
#' @importFrom SeuratObject CreateAssayObject
#' @importFrom Matrix rowSums colSums
#' @concept assay
#'
#' @return A \code{CellGraphAssay} object
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
#' # Create CellGraphAssay
#' cg_assay <- CreateCellGraphAssay(counts = counts, cellgraphs = bipartite_graphs)
#' cg_assay
#'
#' @export
#'
CreateCellGraphAssay <- function (
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

  # Create the Seurat assay object
  counts <- as(counts, "dgCMatrix")
  seurat.assay <- CreateAssayObject(
    counts = counts,
    ...
  )

  # Convert Seurat Assay to a CellGraphAssay
  cg_assay <- as.CellGraphAssay(
    x = seurat.assay,
    cellgraphs = cellgraphs,
    polarization = polarization,
    colocalization = colocalization,
    fs_map = fs_map
  )
  return(cg_assay)
}


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
#' @method CellGraphs MPXAssay
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
#' # CellGraphs getter CellGraphAssay
#' # ---------------------------------
#'
#' # Create CellGraphAssay
#' cg_assay <- CreateCellGraphAssay(counts = counts, cellgraphs = bipartite_graphs)
#' cg_assay
#'
#' # Get cellgraphs from a CellGraphAssay object
#' CellGraphs(cg_assay)
#'
CellGraphs.MPXAssay <- function (
  object,
  ...
) {
  return(slot(object, name = "cellgraphs"))
}


#' @export
#' @method CellGraphs<- MPXAssay
#' @rdname CellGraphs
#' @concept assay
#' @concept cellgraphs
#'
#' @examples
#'
#' # CellGraphs setter CellGraphAssay
#' # ---------------------------------
#'
#' # Set cellgraphs in a CellGraphAssay object
#' CellGraphs(cg_assay) <- cg_assay@cellgraphs
#'
"CellGraphs<-.MPXAssay" <- function (
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


#' @param object A \code{MPXAssay} object, i.e. one of \code{CellGraphAssay} or
#' \code{CellGraphAssay5}
#' @param new.names A character vector with new cell IDs. The length of the vector
#' must be equal to the number of cells in the object and the names must be unique.
#' @param ... Additional arguments (not used)
#'
#' @importFrom SeuratObject RenameCells LayerData
#' @importFrom arrow schema unify_schemas open_dataset arrow_table write_dataset
#'
#' @describeIn MPXAssay-methods Rename cell IDs of a \code{CellGraphAssay} or
#' \code{CellGraphAssay5} object
#' @method RenameCells MPXAssay
#' @concept assay
#' @docType methods
#'
#' @export
RenameCells.MPXAssay <- function (
  object,
  new.names = NULL,
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

  assay <- as(object, Class = ifelse(is(object, "CellGraphAssay"), "Assay", "Assay5"))
  assay <- RenameCells(assay, new.names = new.names, ...)

  # Recreate the CellGraphAssay object
  as_assay_function <- ifelse(is(object, "CellGraphAssay"), as.CellGraphAssay, as.CellGraphAssay5)
  cg_assay_renamed <- as_assay_function(assay)

  # Rename cellgraphs
  cellgraphs <- slot(object, name = "cellgraphs")
  cellgraphs <- set_names(cellgraphs, nm = new.names)
  slot(cg_assay_renamed, name = "cellgraphs") <- cellgraphs

  # Handle polarization slot
  name_conversion <- tibble(new = new.names, component = orig.names)
  polarization <- slot(object, name = "polarization")
  if (length(polarization) > 0) {
    polarization <- polarization %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(cg_assay_renamed, name = "polarization") <- polarization

  # Handle colocalization slot
  colocalization <- slot(object, name = "colocalization")
  if (length(colocalization) > 0) {
    colocalization <- colocalization %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(cg_assay_renamed, name = "colocalization") <- colocalization

  # Handle fs_map
  fs_map <- slot(object, name = "fs_map")
  if (nrow(fs_map) > 0) {
    fs_map$id_map <- fs_map %>% pull(id_map) %>%
      lapply(function(x) {
        x$current_id <- new.names[x$current_id]
        return(x)
      })
  }
  slot(cg_assay_renamed, name = "fs_map") <- fs_map

  return(cg_assay_renamed)
}

#' @inheritParams RenameCells.MPXAssay
#' @describeIn CellGraphAssay-methods Rename cell IDs of a \code{CellGraphAssay} object
#' @concept assay
#' @method RenameCells CellGraphAssay
#' @docType methods
#' @export
#'
RenameCells.CellGraphAssay <- RenameCells.MPXAssay

#' @inheritParams RenameCells.MPXAssay
#' @describeIn CellGraphAssay5-methods Rename cell IDs of a \code{CellGraphAssay5} object
#' @concept assay
#' @method RenameCells CellGraphAssay5
#' @docType methods
#' @export
#'
RenameCells.CellGraphAssay5 <- RenameCells.MPXAssay


#' @param cellgraphs A list of \code{\link{CellGraph}} objects
#' @param polarization A \code{tbl_df} with polarization scores
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param fs_map A \code{tbl_df} with information pxl file paths,
#' sample IDs and component IDs
#'
#' @import rlang
#'
#' @rdname as.CellGraphAssay
#' @method as.CellGraphAssay Assay
#' @concept assay
#'
#' @return A \code{CellGraphAssay} object
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
#' # Create Assay
#' assay <- CreateAssayObject(counts = counts)
#'
#' # Convert Assay to CellGraphAssay
#' cg_assay <- as.CellGraphAssay(assay, cellgraphs = bipartite_graphs)
#' cg_assay
#'
#' @export
#'
as.CellGraphAssay.Assay <- function (
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
      "'cellgraphs' must be a non-empty list with the same number of elements as the number of columns in the Assay" =
        is.list(cellgraphs) &&
        (length(cellgraphs) == ncol(x)))
    stopifnot(
      "'cellgraphs' names must match colnames of the Assay" =
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

  new_assay <- as(object = x, Class = "CellGraphAssay")

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
  from = "Assay",
  to = "CellGraphAssay",
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
        "Class" = "CellGraphAssay"
      ),
      object.list
    )
    return(do.call(what = "new", args = object.list))
  }
)


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


#' @method PolarizationScores MPXAssay
#'
#' @rdname PolarizationScores
#'
#' @export
#'
PolarizationScores.MPXAssay <- function (
  object,
  ...
) {
  slot(object, name = "polarization")
}


#' @method PolarizationScores<- MPXAssay
#'
#' @rdname PolarizationScores
#'
#' @export
#'
"PolarizationScores<-.MPXAssay" <- function (
  object,
  ...,
  value
) {
  # Validate value
  polarization <- .validate_polarization(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "polarization") <- polarization
  return(object)
}


#' @method ColocalizationScores MPXAssay
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
ColocalizationScores.MPXAssay <- function (
  object,
  ...
) {
  slot(object, name = "colocalization")
}


#' @method ColocalizationScores<- MPXAssay
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
"ColocalizationScores<-.MPXAssay" <- function (
  object,
  ...,
  value
) {
  # Validate value
  colocalization <- .validate_colocalization(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "colocalization") <- colocalization
  return(object)
}


#' @method FSMap MPXAssay
#'
#' @rdname FSMap
#'
#' @export
#'
FSMap.MPXAssay <- function (
  object,
  ...
) {
  slot(object, name = "fs_map")
}


#' @method FSMap<- MPXAssay
#'
#' @rdname FSMap
#'
#' @export
#'
"FSMap<-.MPXAssay" <- function (
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


#' MPXAssay Methods
#'
#' Methods for \code{\link{MPXAssay}} objects for generics defined in other
#' packages
#'
#' @param x A \code{\link{MPXAssay}} object
#' @param features Feature names
#' @param cells Cell names
#' @param y A \code{\link{MPXAssay}} object or a list of \code{\link{MPXAssay}} objects
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization);
#' this is recommended if the same normalization approach was applied to all objects
#' @param ... Arguments passed to other methods
#'
#' @name MPXAssay-methods
#' @rdname MPXAssay-methods
#'
#' @concept assay
#'
NULL

#' CellGraphAssay Methods
#'
#' Methods for \code{\link{CellGraphAssay}} objects for generics defined in other
#' packages
#'
#' @param x A \code{\link{CellGraphAssay}} object
#' @param features Feature names
#' @param cells Cell names
#' @param y A \code{\link{CellGraphAssay}} object or a list of \code{\link{CellGraphAssay}} objects
#' @param merge.data Merge the data slots instead of just merging the counts (which requires renormalization);
#' this is recommended if the same normalization approach was applied to all objects
#' @param ... Arguments passed to other methods
#'
#' @name CellGraphAssay-methods
#' @rdname CellGraphAssay-methods
#'
#' @concept assay
#'
NULL

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


#' @describeIn MPXAssay-methods Show method for a \code{CellGraphAssay} or a
#' \code{CellGraphAssay5} object
#' @param object A \code{CellGraphAssay} or a \code{CellGraphAssay5} object
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
#' # Create CellGraphAssay
#' cg_assay <- CreateCellGraphAssay(counts = counts, cellgraphs = bipartite_graphs)
#' cg_assay
#'
#' # Show method
#' cg_assay
#'
#' @export
#'
setMethod (
  f = "show",
  signature = "MPXAssay",
  definition = function(object) {
    cellgraphs <- slot(object, "cellgraphs")
    loaded_graphs <- !sapply(cellgraphs, is.null)
    msg <-
      capture.output(show(as(
        object, Class = ifelse(
          is(object, "CellGraphAssay5"),
          "Assay5",
          "Assay"
        )
      )))
    msg[1] <- gsub(pattern = "^Assay",
                   x = msg[1],
                   replacement = "CellGraphAssay")
    msg <- c(msg, paste0("Loaded CellGraph objects:\n ",
                         sum(loaded_graphs),
                         "\n"))
    cat(paste(msg, collapse = "\n"))
  }
)

#' @describeIn CellGraphAssay-methods Show method for
#' \code{CellGraphAssay} objects
#' @concept assay
#' @method show CellGraphAssay
#' @docType methods
#'
#' @export
#'
setMethod (
  f = "show",
  signature = "CellGraphAssay",
  definition = function(object) as(getMethod("show", "MPXAssay"), "function")(object)
)

#' @describeIn CellGraphAssay5-methods Show method for
#' \code{CellGraphAssay5} objects
#' @concept assay
#' @method show CellGraphAssay5
#' @docType methods
#'
#' @export
#'
setMethod (
  f = "show",
  signature = "CellGraphAssay5",
  definition = function(object) as(getMethod("show", "MPXAssay"), "function")(object)
)


#' @describeIn MPXAssay-methods Subset a \code{CellGraphAssay} or a
#' \code{CellGraphAssay5} object
#' @concept assay
#' @method subset MPXAssay
#' @docType methods
#'
#' @return A \code{MPXAssay} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#' options(Seurat.object.assay.version = "v3")
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' seur <- ReadMPX_Seurat(pxl_file)
#' seur <- LoadCellGraphs(seur)
#' cg_assay <- seur[["mpxCells"]]
#'
#' # Subset CellGraphAssay(5)
#' # ---------------------------------
#' cg_assay_subset <- subset(cg_assay, cells = colnames(cg_assay)[1:3])
#'
#' # Subset Seurat object containing a CellGraphAssay(5)
#' # --------------------------------
#' seur_subset <- subset(seur, cells = colnames(seur)[1:3])
#'
#' @export
#'
subset.MPXAssay <- function (
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
  assay <- as(object = x, Class = ifelse(is(x, "CellGraphAssay"),
                                         "Assay",
                                         "Assay5"))
  assay_subset <- subset(x = assay, features = features, cells = cells)

  # Filter cellgraphs
  cellgraphs_filtered <- cellgraphs[colnames(assay_subset)]

  # Filter cellgraphs
  cellgraphs_filtered <- cellgraphs[colnames(assay_subset)]

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

  # convert standard assay to CellGraphAssay or CellGraphAssay5
  as_assay_func <- ifelse(is(x, "CellGraphAssay"),
                          as.CellGraphAssay,
                          as.CellGraphAssay5)
  cg_assay <- as_assay_func(
    x = assay_subset,
    cellgraphs = cellgraphs_filtered,
    polarization = polarization,
    colocalization = colocalization,
    fs_map = fs_map
  )
  return(cg_assay)
}

#' @describeIn CellGraphAssay-methods Subset a \code{CellGraphAssay} object
#' @concept assay
#' @method subset CellGraphAssay
#' @docType methods
#' @export
#'
subset.CellGraphAssay <- subset.MPXAssay

#' @describeIn CellGraphAssay5-methods Subset a \code{CellGraphAssay5} object
#' @concept assay
#' @method subset CellGraphAssay5
#' @docType methods
#' @export
#'
subset.CellGraphAssay5 <- subset.MPXAssay


#' @param add.cell.ids A character vector with sample names
#'
#' @importFrom SeuratObject Key Key<- RenameCells
#' @describeIn MPXAssay-methods Merge two or more \code{CellGraphAssay} or
#' \code{CellGraphAssay5} objects together
#' @concept assay
#' @method merge MPXAssay
#' @docType methods
#'
#' @importFrom SeuratObject RowMergeSparseMatrices Cells Key Key<- RenameCells
#' @importFrom stringr str_c
#' @importFrom arrow open_dataset write_parquet
#'
#' @examples
#' # Merge multiple CellGraphAssay(5) objects
#' # ---------------------------------
#'
#' # Merge 3 CellGraphAssay(5) objects
#' cg_assay_merged <- merge(cg_assay,
#'                          y = list(cg_assay, cg_assay),
#'                          add.cell.ids = c("A", "B", "C"))
#'
#' @export
#'
merge.MPXAssay <- function (
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  add.cell.ids = NULL,
  ...
) {

  # Validate input parameters
  if (!inherits(y, what = c("list", "MPXAssay")))
    abort("'y' must be a 'CellGraphAssay(5)' object or a list of 'CellGraphAssay(5)' objects")
  if (is.list(y)) {
    for (i in seq_along(y)) {
      if (!inherits(y[[i]], what = "MPXAssay")) {
        abort(glue("Element {i} in 'y' is not a 'CellGraphAssay(5)' object"))
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
  standardassays <- lapply(objects, function(cg_assay) {
    assay <- as(object = cg_assay, Class = ifelse(is(cg_assay, "CellGraphAssay"), "Assay", "Assay5"))
    if (length(Key(assay)) == 0) Key(assay) <- "px_"
    return(assay)
  })

  # Merge Seurat assays
  new_assay <- merge(x = standardassays[[1]],
                     y = standardassays[-1],
                     merge.data = merge.data,
                     ...)

  # Fetch cellgraphs
  cellgraphs_new <- Reduce(c, lapply(objects, function(cg_assay) {
    return(cg_assay@cellgraphs)
  })) %>% set_names(nm = colnames(new_assay))


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
  for (cg_assay in objects) {
    if (length(slot(cg_assay, name = "fs_map")) == 0) {
      fs_map <- NULL
      break
    }
  }
  if (!is.null(fs_map)) {
    fs_map <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
      slot(objects[[i]], name = "fs_map")
    })) %>% mutate(sample = seq_len(n()))
  }

  # Convert to CellGraphAssay and add cellgraphs
  as_assay_func <-
    ifelse(is(cg_assay, "CellGraphAssay"),
           as.CellGraphAssay,
           as.CellGraphAssay5)
  merged_cg_assay <- as_assay_func(
    x = new_assay,
    cellgraphs = cellgraphs_new,
    polarization = polarization,
    colocalization = colocalization,
    fs_map = fs_map
  )

  return(merged_cg_assay)
}

#' @describeIn CellGraphAssay-methods Merge two or more \code{CellGraphAssay} objects together
#' @concept assay
#' @method merge CellGraphAssay
#' @docType methods
#' @export
#'
merge.CellGraphAssay <- merge.MPXAssay

#' @describeIn CellGraphAssay5-methods Merge two or more \code{CellGraphAssay5} objects together
#' @concept assay
#' @method merge CellGraphAssay5
#' @docType methods
#' @export
#'
merge.CellGraphAssay5 <- merge.MPXAssay
