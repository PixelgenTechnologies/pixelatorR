#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod slot slot<- new as slotNames
#' @importClassesFrom Matrix dgCMatrix
NULL

# Declarations used in package check
globalVariables(
  names = c('edges', 'component_new', 'component', 'sample', 'marker_1', 'marker_2'),
  package = 'pixelatorR',
  add = TRUE
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' The CellGraph class
#'
#' The CellGraph class is designed to hold information needed for working with
#' mpx single-cell graphs.
#'
#' @slot cellgraph A \code{tbl_graph} object corresponding to a cell graph
#' @slot counts A \code{matrix}-like object with marker counts
#' @slot layout A \code{list} of \code{tbl_df} objects with coordinates for cell layouts
#'
#' @name CellGraph-class
#' @rdname CellGraph-class
#' @exportClass CellGraph
#' @concept cellgraph
CellGraph <- setClass(
  Class = "CellGraph",
  slots = list(
    cellgraph = "ANY",
    counts = "ANY",
    layout = "ANY"
  )
)


#' The CellGraphAssay class
#'
#' The CellGraphAssay object is an extended \code{\link[SeuratObject]{Assay}}
#' for the storage and analysis of mpx single-cell data.
#'
#' @slot cellgraphs A named list of \code{\link{CellGraph}} objects
#' @slot polarization A \code{tbl_df} with polarization scores
#' @slot colocalization A \code{tbl_df} with colocalization scores
#' @slot arrow_dir A character giving the name of the directory where the edgelist
#' parquet file(s) are stored
#' @slot arrow_data An \code{R6} class object with an arrow Dataset
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
    arrow_dir = "character",
    arrow_data = "ANY"
  )
)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a CellGraph object
#'
#' @param cellgraph A \code{tbl_graph} object representing an mpx single-cell graph
#' @param counts A \code{dgCMatrix} with marker counts
#' @param layout A \code{tbl_df} object with cell layout(s)
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @concept cellgraph
#'
#' @return A \code{CellGraph} object
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' edge_list <-
#' ReadMPX_item(
#'   system.file("extdata/PBMC_10_cells",
#'               "Sample01_test.pxl",
#'               package = "pixelatorR"),
#'   items = "edgelist"
#' )
#' bipart_graph <-
#'   edge_list %>%
#'   select(upia, upib, marker) %>%
#'   distinct() %>%
#'   as_tbl_graph(directed = FALSE) %>%
#'   mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
#' attr(bipart_graph, "type") <- "bipartite"
#'
#' cg <- CreateCellGraphObject(cellgraph = bipart_graph)
#' cg
#'
#' @export
#'
CreateCellGraphObject <- function (
  cellgraph,
  counts = NULL,
  layout = NULL,
  verbose = FALSE
) {

  # Check input parameters
  stopifnot("'cellgraph' must be a non-empty 'tbl_graph' object" = inherits(cellgraph, what = "tbl_graph") & (length(cellgraph) > 0))
  if (!is.null(counts)) {
    stopifnot("'counts' must be a non-empty 'dgCMatrix' object" = inherits(counts, what = "dgCMatrix") & (length(counts) > 0))
  }
  if (!is.null(layout)) {
    stopifnot("'layout' must be a non-empty 'tbl_df' object" = inherits(layout, what = "tbl_df") & (length(layout) > 0))
  }

  if (!"type" %in% names(attributes(cellgraph))) {
    abort("Graph attribute 'type' is missing.")
  } else {
    if (verbose && check_global_verbosity())
      cli_alert_info("Got a graph of type '{attr(cellgraph, 'type')}'")
  }

  # Add checks for graph types
  if (attr(cellgraph, "type") == "bipartite") {
    stopifnot("Node attribute 'name' is missing" = "name" %in% vertex_attr_names(cellgraph))
    stopifnot("Node attribute 'node_type' is missing" = "node_type" %in% vertex_attr_names(cellgraph))
  }
  #TODO: Add check for A-node-projection and linegraph

  # create object
  cellgraph <- new(
    Class = "CellGraph",
    cellgraph = cellgraph,
    counts = counts,
    layout = layout
  )

  return(cellgraph)
}

#' Create a CellGraphAssay object
#'
#' Create a \code{\link{CellGraphAssay}} object from a count matrix. The expected
#' format of the input matrix is features x cells.
#'
#' @param counts Unnormalized data (raw counts)
#' @param cellgraphs A named list of \code{\link{CellGraph}} objects
#' @param polarization A \code{tbl_df} with polarization scores
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param arrow_dir A path to an existing directory
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
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
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
  arrow_dir = NULL,
  outdir = NULL,
  overwrite = FALSE,
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

  # Set path to NA if NULL
  arrow_dir <- arrow_dir %||% NA_character_

  # Validate polarization and colocalization
  polarization <- .validate_polarization(polarization, cell_ids = colnames(counts), markers = rownames(counts), verbose = verbose)
  colocalization <- .validate_colocalization(colocalization, cell_ids = colnames(counts), markers = rownames(counts), verbose = verbose)

  # If arrow_dir are provided, attempt to lazy load edgelist
  if (!is.na(arrow_dir)) {
    stopifnot("'arrow_dir' must be a non-empty character" = is.character(arrow_dir) && (length(arrow_dir) >= 1))
    for (path in arrow_dir) {
      if (!(dir.exists(path) || file.exists(path))) {
        abort(glue("Directory/file {path} doesn't exist"))
      }
    }
    arrow_dir <- sapply(arrow_dir, normalizePath)

    # Load arrow dataset
    fsd <- ReadMPX_arrow_edgelist(path = arrow_dir, outdir = outdir, return_list = TRUE, overwrite = overwrite, verbose = FALSE)

    # Update arrow_dir to new temporary directory
    arrow_dir <- fsd$arrow_dir
    fsd <- fsd$ArrowObject
    if (!all(c("upia", "marker", "component", "sample") %in% names(fsd))) {
      abort(glue("The following columns need to be present in the edgelist:",
                 " 'upia', 'marker', 'component', 'sample'\n",
                 "Check that {arrow_dir} is correctly formatted."))
    }

    # check for required fields
    stopifnot("column 'component' is missing from cellgraphs" = "component" %in% names(fsd))
    available_components <- fsd %>% pull(component, as_vector = TRUE) %>% unique()
    if (!all(colnames(counts) %in% available_components))
      abort(glue("Some components are not available in the edge list"))
  } else {
    fsd <- NULL
  }

  # Create the Seurat assay
  seurat.assay <- CreateAssayObject(
    counts = counts,
    ...
  )

  # Convert Seurat Assay to a CellGraphAssay
  pxCell.assay <- as.CellGraphAssay(
    x = seurat.assay,
    cellgraphs = cellgraphs,
    polarization = polarization,
    colocalization = colocalization,
    arrow_dir = arrow_dir,
    arrow_data = fsd
  )
  return(pxCell.assay)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Get methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get and set CellGraph object data
#'
#' @param object A \code{\link{CellGraph}} object
#' @param slot Information to pull from object (cellgraph, meta_data, layout)
#'
#' @import rlang
#'
#' @rdname CellGraphData
#'
#' @return \code{GetCellGraphData}: A \code{\link{CellGraph}} object slot
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' edge_list <-
#'   ReadMPX_item(
#'     system.file("extdata/PBMC_10_cells",
#'                 "Sample01_test.pxl",
#'                 package = "pixelatorR"),
#'     items = "edgelist"
#'   )
#' bipart_graph <-
#'   edge_list %>%
#'   select(upia, upib, marker) %>%
#'   distinct() %>%
#'   as_tbl_graph(directed = FALSE) %>%
#'   mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
#' attr(bipart_graph, "type") <- "bipartite"
#'
#' cg <- CreateCellGraphObject(cellgraph = bipart_graph)
#'
#' # Get slot data
#' CellGraphData(cg, slot = "cellgraph")
#'
#' @export
#'
CellGraphData <- function (
  object,
  slot = "cellgraph"
) {
  if (!inherits(object, what = "CellGraph")) abort(glue("Invalid class {class(object)}"))
  if (!(slot %in% slotNames(x = object))) {
    abort(glue("slot must be one of {paste(slotNames(x = object), collapse = ', ')}"))
  }
  return(slot(object = object, name = slot))
}


#' @param value A new variable to place in \code{slot}
#'
#' @rdname CellGraphData
#'
#' @return \code{CellGraphData<-}: A \code{\link{CellGraph}} with updated data
#'
#' @examples
#' # Set slot data
#' CellGraphData(cg, slot = "cellgraph") <- CellGraphData(cg, slot = "cellgraph")
#'
#'
#' @export
#'
"CellGraphData<-" <- function (
  object,
  slot = "cellgraph",
  value
) {
  if (!inherits(object, what = "CellGraph")) abort(glue("Invalid class {class(object)}"))
  if (!(slot %in% slotNames(x = object))) {
    abort(glue("slot must be one of {paste(slotNames(x = object), collapse = ', ')}"))
  }

  # Get counts and layouts
  cellgraph <- slot(object, name = "cellgraph")
  counts <- slot(object, name = "counts")
  layouts <- slot(object, name = "layout")

  # Validate cellgraph input
  if (slot == "cellgraph") {
    stopifnot(
      "'value' must be a 'tbl_graph' object" =
        inherits(value, what = "tbl_graph")
    )
    if (length(counts) > 0) {
      stopifnot(
        "names in 'value' do not match the row names of 'counts'" =
          all(names(value) == rownames(counts))
      )
    }
    if (length(layouts) > 0) {
      for (layout in names(layouts)) {
        if (length(value) != nrow(layouts[[layout]])) {
          abort(glue("Number of nodes in 'value' does not match the ",
                     "number of rows in the '{layout}' layout table"))
        }
      }
    }
    slot(object, name = "cellgraph") <- cellgraph
  }

  # Validate counts input
  if (slot == "counts") {
    stopifnot(
      "'value' must be a 'dgCMatrix' object" =
        inherits(value, what = "dgCMatrix")
    )
    if (length(layouts) > 0) {
      for (layout in names(layouts)) {
        if (nrow(counts) != nrow(layouts[[layout]])) {
          abort(glue("Number of rows in 'value' does not match the ",
                     "number of rows in the '{layout}' layout table"))
        }
      }
    }
    stopifnot(
      "Number of rows in 'value' do not match the number of nodes in 'cellgraph'" =
        length(cellgraph) == nrow(value)
    )
    slot(object, name = "counts") <- counts
  }

  # Validate layout input
  if (slot == "layout") {
    stopifnot(
      "'value' must be a 'list'" =
        inherits(value, what = "list"),
      "'value' is empty" =
        length(value) > 0,
      "'names' of 'value' cannot be NULL" =
        !is.null(names(value))
    )
    for (layout in names(value)) {
      stopifnot(
        inherits(value[[layout]], what = "tbl_df")
      )
      if (length(cellgraph) != nrow(value[[layout]])) {
        abort(glue("Number of nodes in 'cellgraph' does not match the ",
                   "number of rows in the '{layout}' layout table provided in 'value'"))
      }
      if (length(counts) > 0) {
        if (nrow(counts) != nrow(value[[layout]])) {
          abort(glue("Number of rows in 'counts' does not match ",
                     "the number of rows in '{layout}' table provided in 'value'"))
        }
      }
    }
    slot(object, name = "layout") <- value
  }

  return(object)
}

#' @rdname CellGraphs
#' @method CellGraphs CellGraphAssay
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
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
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
CellGraphs.CellGraphAssay <- function (
  object,
  ...
) {
  return(slot(object, name = "cellgraphs"))
}

#' @param object A Seurat object or CellGraphAssay object
#' @importFrom SeuratObject DefaultAssay
#' @rdname CellGraphs
#' @method CellGraphs Seurat
#' @export
#' @concept assay
#' @concept cellgraphs
#'
#' @examples
#'
#' # CellGraphs getter Seurat
#' # ---------------------------------
#' se <- ReadMPX_Seurat(pxl_file, overwrite = TRUE, return_cellgraphassay = TRUE)
#'
#' # Get cellgraphs from a Seurat object
#' CellGraphs(se)
#'
#'
CellGraphs.Seurat <- function (
  object,
  ...
) {
  assay <- DefaultAssay(object = object)
  return(CellGraphs(object = object[[assay]]))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @export
#' @method CellGraphs<- CellGraphAssay
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
"CellGraphs<-.CellGraphAssay" <- function (
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

#' @export
#' @method CellGraphs<- Seurat
#' @rdname CellGraphs
#' @concept assay
#' @concept cellgraphs
#' @importFrom SeuratObject DefaultAssay
#'
#' @examples
#' # CellGraphs setter Seurat
#' # ---------------------------------
#'
#' # Set cellgraphs in a Seurat object
#' CellGraphs(se) <- cg_assay@cellgraphs
#'
"CellGraphs<-.Seurat" <- function (
  object,
  ...,
  value
) {
  assay <- DefaultAssay(object = object)
  CellGraphs(object = object[[assay]]) <- value
  return(object)
}

#' @importFrom SeuratObject RenameCells LayerData
#' @importFrom arrow schema unify_schemas open_dataset arrow_table write_dataset
#'
#' @rdname RenameCells
#' @method RenameCells CellGraphAssay
#' @concept assay
#'
#' @export
RenameCells.CellGraphAssay <- function (
  object,
  new.names = NULL,
  ...
) {

  stopifnot(
    "'new.names' must be a character vector with the same length as the number of cells present in 'object'" =
      inherits(new.names, what = "character") & (length(new.names) == ncol(object))
  )
  names(slot(object = object, name = "cellgraphs")) <- new.names

  # save original cell IDs
  orig.names <- colnames(object)

  names(x = new.names) <- NULL
  for (data.slot in object[]) {
    old.data <- LayerData(object = object, layer = data.slot)
    if (ncol(x = old.data) <= 1) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }

  # Handle arrow dataset
  if (!is.na(slot(object, name = "arrow_dir"))) {
    # Restore arrow connection if broken
    object <- RestoreArrowConnection(object, verbose = FALSE)

    # Get schema for the component column (string or large string)
    schema <- unify_schemas(schema(object@arrow_data)["component"],
                            schema(object@arrow_data)["component"] %>% setNames(nm = "component_new"))

    # Create conversion table with defined schema
    conv_table <- arrow_table(component = orig.names, component_new = new.names, schema = schema)

    fsd <- slot(object, name = "arrow_data")

    # Left join new names and save edgelist
    session_tmpdir_random <- file.path(getOption("pixelatorR.arrow_outdir"), paste0(.generate_random_string(), "-", format(Sys.time(), "%Y-%m-%d-%H%M%S")))

    # arrow doesn't support pasting to create new character columns,
    # instead we do a left join with our conversion table to
    # obtain the new cell/component names
    fsd %>%
      left_join(y = conv_table, by = "component") %>% # Add new cell names
      select(-component) %>% # remove column with old names
      rename(component = component_new) %>% # rename column
      group_by(sample) %>%
      write_dataset(path = session_tmpdir_random) # write data to a parquet file

    # Rename parquet file to ensure consistency with other functions
    files <- list.files(session_tmpdir_random, pattern = "parquet", recursive = TRUE, full.names = TRUE)
    for (f in files) {
      if (basename(f) != "edgelist.parquet") {
        file.rename(from = f, to = file.path(dirname(f), "edgelist.parquet"))
      }
    }

    # Reload arrow dataset
    fsd_new <- open_dataset(session_tmpdir_random)
    slot(object, name = "arrow_data") <- fsd_new
    slot(object, name = "arrow_dir") <- session_tmpdir_random
  }

  # Handle polarization slot
  name_conversion <- tibble(new = new.names, component = orig.names)
  polarization <- slot(object, name = "polarization")
  if (length(polarization) > 0) {
    polarization <- polarization %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(object, name = "polarization") <- polarization

  # Handle colocalization slot
  colocalization <- slot(object, name = "colocalization")
  if (length(colocalization) > 0) {
    colocalization <- colocalization %>%
      left_join(y = name_conversion, by = "component") %>%
      select(-component) %>%
      rename(component = new)
  }
  slot(object, name = "colocalization") <- colocalization

  return(object)
}


#' @param cellgraphs A list of \code{\link{CellGraph}} objects
#' @param polarization A \code{tbl_df} with polarization scores
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param arrow_dir TODO
#' @param arrow_data An \code{R6} class object created with arrow
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
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
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
  arrow_dir = NULL,
  arrow_data = NULL,
  ...
) {

  # Check cellgraphs
  if (!is.null(cellgraphs)) {
    stopifnot("'cellgraphs' must be a non-empty list with the same number of elements as the number of columns in the Assay" = is.list(cellgraphs) & (length(cellgraphs) == ncol(x)))
    stopifnot("'cellgraphs' names must match colnames of the Assay" = all(names(cellgraphs) == colnames(x)))
    for (i in seq_along(cellgraphs)) {
      if (!inherits(x = cellgraphs[[i]], what = c("CellGraph", "NULL"))) {
        abort(glue("Element {i} is not a CellGraph object or NULL"))
      }
    }
  } else {
    cellgraphs <- rep(list(NULL), ncol(x)) %>% setNames(nm = colnames(x))
  }

  # Validate polarization and colocalization
  polarization <- .validate_polarization(polarization, cell_ids = colnames(x), markers = rownames(x))
  colocalization <- .validate_colocalization(colocalization, cell_ids = colnames(x), markers = rownames(x))

  # Abort if cellgraphs is empty and neither arrow_dir or arrow_data is provided
  loaded_graphs <- sum(sapply(cellgraphs, is.null))
  if (loaded_graphs == ncol(x)) {
    stopifnot("One of 'arrow_dir' or 'arrow_data' must be provided if 'cellgraphs is empty'" = (!is.null(arrow_dir)) || (!is.null(arrow_data)))
  }

  # Handle arrow_dir
  arrow_dir <- arrow_dir %||% NA_character_
  if (!is.na(arrow_dir)) {
    stopifnot("'arrow_dir' must be a non-empty character" = is.character(arrow_dir) && (length(arrow_dir) >= 1))
    for (path in arrow_dir) {
      if (!(dir.exists(path) || file.exists(path))) {
        abort(glue("Directory/file {path} doesn't exist"))
      }
    }
    arrow_dir <- sapply(arrow_dir, normalizePath) %>% unname()
  }

  # Handle arrow_data
  if (!is.null(arrow_data)) {
    #TODO: validate class
    stopifnot(all(c("upia", "marker", "component", "sample") %in% names(arrow_data)))
    stopifnot("A valid 'arrow_dir' must be provided if 'arrow_data' is provided" = !is.na(arrow_dir))
    stopifnot("column 'component' is missing from 'arrow_data" = "component" %in% names(arrow_data))
    stopifnot("One or several components are missing from 'arrow_data'" = all(colnames(x) %in% (arrow_data %>% pull(component, as_vector = TRUE))))
  }

  new.assay <- as(object = x, Class = "CellGraphAssay")
  slot(new.assay, name = "arrow_dir") <- NA_character_

  # Add slots
  slot(new.assay, name = "cellgraphs") <- cellgraphs
  slot(new.assay, name = "polarization") <- polarization
  slot(new.assay, name = "colocalization") <- colocalization
  if (!is.null(arrow_data)) {
    # If an R6 class object exists, add it and arrow_dir to the appropriate slots
    slot(new.assay, name = "arrow_data") <- arrow_data
    slot(new.assay, name = "arrow_dir") <- arrow_dir
  }

  return(new.assay)
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

#' @method PolarizationScores CellGraphAssay
#'
#' @rdname PolarizationScores
#'
#' @export
#'
PolarizationScores.CellGraphAssay <- function (
  object,
  ...
) {
  slot(object, name = "polarization")
}

#' @param assay Name of a \code{CellGraphAssay}
#' @param meta_data_columns A character vector with meta.data column names.
#' This option can be useful to join meta.data columns with the polarization
#' score table.
#'
#' @method PolarizationScores Seurat
#'
#' @rdname PolarizationScores
#'
#' @export
#'
PolarizationScores.Seurat <- function (
  object,
  assay = NULL,
  meta_data_columns = NULL,
  ...
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("Assay '{assay}' is not a CellGraphAssay"))
  }

  # Get polarizaation scores from CellGraphAssay
  pol_scores <- PolarizationScores(object[[assay]])

  # Handle adding meta data columns
  if (!is.null(meta_data_columns)) {
    stopifnot(
      "'meta_data_columns' must be a non-empty character vector" =
        is.character(meta_data_columns) &&
        (length(meta_data_columns) > 0)
    )
    meta_data_columns_valid <- meta_data_columns %in% colnames(object[[]])
    if (any(!meta_data_columns_valid)) {
      abort(glue("The following columns were not found in the meta.data slot: ",
      "{paste(meta_data_columns[!meta_data_columns_valid], collapse=', ')}"))
    }

    # Add additional meta.data slots
    pol_scores <- pol_scores %>%
      left_join(y = object[[]] %>%
                  as_tibble(rownames = "component") %>%
                  select(component, all_of(meta_data_columns)),
                by = "component")
  }

  return(pol_scores)
}

#' @method PolarizationScores<- CellGraphAssay
#'
#' @rdname PolarizationScores
#'
#' @export
#'
"PolarizationScores<-.CellGraphAssay" <- function (
  object,
  ...,
  value
) {
  # Validate value
  polarization <- .validate_polarization(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "polarization") <- polarization
  return(object)
}

#' @method PolarizationScores<- Seurat
#'
#' @rdname PolarizationScores
#'
#' @export
#'
"PolarizationScores<-.Seurat" <- function (
  object,
  assay = NULL,
  ...,
  value
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  # Fetch polarization scores from assay
  cg_assay <- object[[assay]]
  PolarizationScores(cg_assay) <- value
  object[[assay]] <- cg_assay

  return(object)
}

#' @method ColocalizationScores CellGraphAssay
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
ColocalizationScores.CellGraphAssay <- function (
  object,
  ...
) {
  slot(object, name = "colocalization")
}

#' @param assay Name of a \code{CellGraphAssay}
#' @param meta_data_columns A character vector with meta.data column names.
#' This option can be useful to join meta.data columns with the polarization
#' score table.
#'
#' @method ColocalizationScores Seurat
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
ColocalizationScores.Seurat <- function (
  object,
  assay = NULL,
  meta_data_columns = NULL,
  ...
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  cg_assay <- object[[assay]]
  if (!inherits(cg_assay, what = "CellGraphAssay")) {
    abort(glue("Assay '{assay}' is not a CellGraphAssay"))
  }

  # Get colocalization scores
  coloc_scores <- ColocalizationScores(cg_assay)

  # Handle adding meta data columns
  if (!is.null(meta_data_columns)) {
    stopifnot(
      "'meta_data_columns' must be a non-empty character vector" =
        is.character(meta_data_columns) &&
        (length(meta_data_columns) > 0)
    )
    meta_data_columns_valid <- meta_data_columns %in% colnames(object[[]])
    if (any(!meta_data_columns_valid)) {
      abort(glue("The following columns were not found in the meta.data slot: ",
                 "{paste(meta_data_columns[!meta_data_columns_valid], collapse=', ')}"))
    }

    # Add additional meta.data slots
    coloc_scores <- coloc_scores %>%
      left_join(y = object[[]] %>%
                  as_tibble(rownames = "component") %>%
                  select(component, all_of(meta_data_columns)),
                by = "component")
  }

  return(coloc_scores)
}

#' @method ColocalizationScores<- CellGraphAssay
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
"ColocalizationScores<-.CellGraphAssay" <- function (
  object,
  ...,
  value
) {
  # Validate value
  colocalization <- .validate_colocalization(value, cell_ids = colnames(object), markers = rownames(object))
  slot(object, name = "colocalization") <- colocalization
  return(object)
}

#' @method ColocalizationScores<- Seurat
#'
#' @rdname ColocalizationScores
#'
#' @export
#'
"ColocalizationScores<-.Seurat" <- function (
  object,
  assay = NULL,
  ...,
  value
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  # Fetch polarization scores from assay
  cg_assay <- object[[assay]]
  ColocalizationScores(cg_assay) <- value
  object[[assay]] <- cg_assay

  return(object)
}


#' @method ArrowData CellGraphAssay
#'
#' @rdname ArrowData
#'
#' @export
#'
ArrowData.CellGraphAssay <- function (
  object,
  ...
) {
  slot(object, name = "arrow_data")
}


#' @param assay Name of a \code{CellGraphAssay}
#'
#' @method ArrowData Seurat
#'
#' @rdname ArrowData
#'
#' @export
#'
ArrowData.Seurat <- function (
  object,
  assay = NULL,
  ...
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  ArrowData(object[[assay]])
}


#' @method ArrowData<- CellGraphAssay
#'
#' @rdname ArrowData
#'
#' @export
#'
"ArrowData<-.CellGraphAssay" <- function (
  object,
  ...,
  value
) {
  # Validate value
  if (!inherits(value, what = c("Dataset", "arrow_dplyr_query", "ArrowObject"))) {
    abort(glue("Invalid class '{class(value)[1]}'."))
  }
  msg <- tryCatch(value %>% nrow(), error = function(e) "Error")
  if (msg == "Error") {
    abort("Invalid 'value'")
  }
  slot(object = object, name = "arrow_data") <- value
  return(object)
}


#' @method ArrowData<- Seurat
#'
#' @rdname ArrowData
#'
#' @export
#'
"ArrowData<-.Seurat" <- function (
  object,
  assay = NULL,
  ...,
  value
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  ArrowData(cg_assay) <- value
  object[[assay]] <- cg_assay

  return(object)
}


#' @param assay Name of a \code{CellGraphAssay}
#'
#' @method ArrowDir CellGraphAssay
#'
#' @rdname ArrowDir
#'
#' @export
#'
ArrowDir.CellGraphAssay <- function (
  object,
  ...
) {
  slot(object, name = "arrow_dir")
}

#' @method ArrowDir Seurat
#'
#' @rdname ArrowDir
#'
#' @export
#'
ArrowDir.Seurat <- function (
  object,
  assay = NULL,
  ...
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)

  return(ArrowDir(object[[assay]]))
}


#' @method ArrowDir<- CellGraphAssay
#'
#' @rdname ArrowDir
#'
#' @export
#'
"ArrowDir<-.CellGraphAssay" <- function (
  object,
  ...,
  value
) {
  stopifnot(
    "'value' must be a non-empty character vector" =
      is.character(value) &&
      (length(value) == 1)
  )
  if (!file.exists(value)) {
    abort(glue("{value} doesn't exist"))
  }
  slot(object = object, name = "arrow_dir") <- value
  return(object)
}


#' @method ArrowDir<- Seurat
#'
#' @rdname ArrowDir
#'
#' @export
#'
"ArrowDir<-.Seurat" <- function (
  object,
  assay = NULL,
  ...,
  value
) {

  # Use default assay if assay = NULL
  assay <- assay %||% DefaultAssay(object)
  cg_assay <- object[[assay]]
  ArrowDir(cg_assay) <- value
  object[[assay]] <- cg_assay

  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for R-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Base methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Show method for \code{CellGraph} object
#'
#' @param object A \code{CellGraph} object
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' edge_list <-
#' ReadMPX_item(
#'   system.file("extdata/PBMC_10_cells",
#'               "Sample01_test.pxl",
#'               package = "pixelatorR"),
#'   items = "edgelist"
#' )
#' bipart_graph <-
#'   edge_list %>%
#'   select(upia, upib, marker) %>%
#'   distinct() %>%
#'   as_tbl_graph(directed = FALSE) %>%
#'   mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
#' attr(bipart_graph, "type") <- "bipartite"
#'
#' cg <- CreateCellGraphObject(cellgraph = bipart_graph)
#'
#' # Show method
#' cg
#'
setMethod (
 f = "show",
 signature = "CellGraph",
 definition = function(object) {
   graph_type <- attr(slot(object, "cellgraph"), "type")
   if (is.null(slot(object, "counts"))) {
     n_markers <- NULL
   } else {
     n_markers <- ncol(slot(object, "counts"))
   }
   cat(
     "A CellGraph object containing a", graph_type, "graph with",
     slot(object = object, name = "cellgraph") %>% length() %>% col_br_blue(),
     "nodes and",
     slot(object = object, name = "cellgraph") %>% gsize() %>% col_br_blue(),
     "edges "
   )
   if (is.null(n_markers)) {
     cat("\n")
   } else {
     cat("covering", n_markers, "markers\n")
   }
 }
)

#' Show method for \code{CellGraphAssay} object
#'
#' @param object A \code{CellGraphAssay} object
#'
#' @importFrom methods show
#'
#' @examples
#'
#' library(pixelatorR)
#' library(dplyr)
#' library(tidygraph)
#'
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
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
setMethod (
  f = "show",
  signature = "CellGraphAssay",
  definition = function(object) {
    cellgraphs <- slot(object, "cellgraphs")
    loaded_graphs <- !sapply(cellgraphs, is.null)
    show(as(object, Class = "Assay"))
    cat(
      "Loaded CellGraph objects:",
      sum(loaded_graphs),
      "\n"
    )
  }
)

#' @describeIn CellGraphAssay-methods Subset a \code{CellGraphAssay} object
#' @importClassesFrom SeuratObject Assay
#' @concept assay
#' @method subset CellGraphAssay
#'
#' @importFrom arrow write_dataset
#'
#' @return A \code{CellGraphAssay} object
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' pxl_file <- system.file("extdata/PBMC_10_cells",
#'                         "Sample01_test.pxl",
#'                         package = "pixelatorR")
#' seur <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
#' seur <- LoadCellGraphs(seur)
#' cg_assay <- seur[["mpxCells"]]
#'
#' # Subset CellGraphAssay
#' # ---------------------------------
#' cg_assay_subset <- subset(cg_assay, cells = colnames(cg_assay)[1:3])
#'
#' # Subset Seurat object containing a CellGraphAssay
#' # --------------------------------
#' seur_subset <- subset(seur, cells = colnames(cg_assay)[1:3])
#'
#' # Compare size of edge lists stored on disk
#' ArrowData(seur) %>% dim()
#' ArrowData(seur_subset) %>% dim()
#'
#' @export
#'
subset.CellGraphAssay <- function (
  x,
  features = NULL,
  cells = NULL,
  ...
) {

  # Get cellgraphs
  cellgraphs <- x@cellgraphs

  # subset elements in the standard assay
  standardassay <- as(object = x, Class = "Assay")
  standardassay <- subset(x = standardassay, features = features, cells = cells)

  # Filter cellgraphs
  cellgraphs_filtered <- cellgraphs[colnames(standardassay)]

  # Fetch arrow_dir
  arrow_dir <- slot(x, name = "arrow_dir")

  # Handle arrow data set if available
  if (!is.na(arrow_dir)) {
    x <- RestoreArrowConnection(x, verbose = FALSE)

    # Filter edgelists
    if (!is.null(cells)) {

      # Only filter by cells
      if (length(cells) < ncol(x)) {

        # Create a temporary directory with a unique name
        session_tmpdir_random <- file.path(getOption("pixelatorR.arrow_outdir"), paste0(.generate_random_string(), "-", format(Sys.time(), "%Y-%m-%d-%H%M%S")))

        # Filter edgelist and export it
        slot(x, name = "arrow_data") %>%
          filter(component %in% cells) %>%
          group_by(sample) %>%
          write_dataset(session_tmpdir_random)

        # Rename parquet files for consistensy with other functions
        files <- list.files(session_tmpdir_random, pattern = "parquet", recursive = TRUE, full.names = TRUE)
        for (f in files) {
          if (basename(f) != "edgelist.parquet") {
            file.rename(from = f, file.path(dirname(f), "edgelist.parquet"))
          }
        }

        # Update arrow_dir
        arrow_dir <- session_tmpdir_random
      }
    }
    slot(x, name = "arrow_data") <- arrow::open_dataset(arrow_dir)
  }

  # Filter polarization and colocaliation
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

  # convert standard assay to CellGraphAssay
  pxcellassay <- as.CellGraphAssay(
    x = standardassay,
    cellgraphs = cellgraphs_filtered,
    polarization = polarization,
    colocalization = colocalization,
    arrow_dir = arrow_dir,
    arrow_data = slot(x, name = "arrow_data")
  )
  return(pxcellassay)
}

#' @describeIn CellGraphAssay-methods Merge two or more \code{CellGraphAssay} objects together
#' @concept assay
#' @method merge CellGraphAssay
#'
#' @importFrom SeuratObject RowMergeSparseMatrices Cells Key Key<- RenameCells
#' @importFrom stringr str_c
#' @importFrom arrow open_dataset write_parquet
#'
#' @examples
#' # Merge CellGraphAssays
#' # ---------------------------------
#'
#' # Merge 3 CellGraphAssays
#' cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay))
#'
#' # Check size of merged edge lists stored on disk
#' ArrowData(cg_assay_merged) %>% dim()
#'
#' @export
#'
merge.CellGraphAssay <- function (
  x = NULL,
  y = NULL,
  merge.data = TRUE,
  ...
) {

  # Validate input parameters
  if (!inherits(y, what = c("list", "CellGraphAssay")))
    abort("'y' must be a 'CellGraphAssay' object or a list of 'CellGraphAssay' objects")
  if (is.list(y)) {
    for (i in seq_along(y)) {
      if (!inherits(y[[i]], what = "CellGraphAssay")) {
        abort(glue("Element {i} in 'y' is not a 'CellGraphAssay'"))
      }
    }
  }

  objects <- c(x, y)

  # Check duplicate cell names
  cell.names <- unlist(lapply(objects, colnames))
  if (any(duplicated(x = cell.names))) {
    cli_alert_warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    objects <-
      lapply(seq_along(objects), function(i) {
        return(RenameCells(object = objects[[i]], new.names = paste0(Cells(x = objects[[i]]), "_", i)))
      })
  }

  # Fetch cellgraphs
  cellgraphs_new <- Reduce(c, lapply(seq_along(objects), function(i) {
    el <- objects[[i]]
    return(el@cellgraphs)
  }))

  # Fetch standard assays
  standardassays <- lapply(objects, function(el) {
    assay <- as(object = el, Class = "Assay")
    if (length(Key(assay)) == 0) Key(assay) <- "px_"
    return(assay)
  })

  # Fetch all arrow_dirs
  all_arrow_dirs <- sapply(objects, function(el) slot(el, name = "arrow_dir"))

  # Merge Seurat assays
  new_assay <- merge(x = standardassays[[1]],
                     y = standardassays[2:length(standardassays)],
                     merge.data = merge.data,
                     ...)

  # Rename cellgraphs
  names(cellgraphs_new) <- colnames(new_assay)

  # Run edgelist merge if arrow_dirs exists
  if (!any(is.na(all_arrow_dirs))) {

    if (!any(dir.exists(all_arrow_dirs)))
      abort(glue("Paths to edgelist parquet file directories {paste(all_arrow_dirs, collapse = ',')} doesn't exist"))

    # Create a new temporary directory with a time stamp
    new_dir <- file.path(getOption("pixelatorR.arrow_outdir"), paste0(.generate_random_string(), "-", format(Sys.time(), "%Y-%m-%d-%H%M%S")))
    dir.create(path = new_dir, showWarnings = FALSE)

    # List hive-partitioned directories for all objects
    dirs_unpacked <- do.call(bind_rows, lapply(seq_along(all_arrow_dirs), function(i) {
      cur_dir <- all_arrow_dirs[i]
      dirs_cur <- list.files(cur_dir)
      cur_unpacked <- do.call(bind_rows, strsplit(dirs_cur, "=") %>%
                                lapply(function(x) {tibble(id = x[1], sample = x[2])}))
      cur_unpacked <- cur_unpacked %>%
        mutate(old_dir = list.files(cur_dir, full.names = TRUE), ID = i)
    }))

    # Create new names
    dirs_unpacked <- dirs_unpacked %>%
      mutate(sample = paste0("S", 1:n())) %>%
      mutate(new_dir = file.path(new_dir, paste0("sample=", sample)))

    # Move parquet files
    for (i in 1:nrow(dirs_unpacked)) {
      unlink(dirs_unpacked$new_dir[i], recursive = TRUE)
      dir.create(path = dirs_unpacked$new_dir[i], showWarnings = TRUE)
      file.copy(from = list.files(dirs_unpacked$old_dir[i], full.names = TRUE, recursive = TRUE),
                to = dirs_unpacked$new_dir[i], recursive = TRUE)
    }

    # Load parquet files
    arrow_data <- open_dataset(new_dir)

    # Change ids in parquet edgelist if necessary
    # This check could become slow for larger datasets
    # TODO: Look for alternative solution
    if (!all(colnames(new_assay) %in% (arrow_data %>% pull(component, as_vector = TRUE)))) {

      # List all parquet files in new_dir
      all_parquet_files <- list.files(new_dir, recursive = TRUE, full.names = TRUE)

      # Update cell IDs with the same rules as Seurat:::CheckDuplicateCellNames
      for (i in seq_along(all_parquet_files)) {
        parquet_file <- all_parquet_files[i]
        edgelist <- open_dataset(parquet_file) %>%
          mutate(id = i) %>%
          mutate(component = str_c(component, id, sep = "_")) %>%
          select(-id)
        write_parquet(x = edgelist, sink = parquet_file)
        rm(edgelist)
      }
      # Load parquet files from new directory
      arrow_data <- open_dataset(new_dir)
    }
  } else {
    # If any path in the input objects is NA, simply inactivate the
    # slots required to handle the arrow Dataset
    arrow_data <- NULL
    new_dir <- NA_character_
  }

  # Merge polarization and colocalization scores
  name_conversion <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
    tibble(new = colnames(objects[[i]]), sample = i)
  })) %>% mutate(component = cell.names) %>%
    group_by(sample) %>%
    group_split()
  polarization <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
    pl <- slot(objects[[i]], name = "polarization")
    if (length(pl) == 0) return(pl)
    pl <- pl %>% left_join(name_conversion[[i]], by = "component") %>%
      select(-component, -sample) %>%
      rename(component = new)
    return(pl)
  }))
  colocalization <- do.call(bind_rows, lapply(seq_along(objects), function(i) {
    cl <- slot(objects[[i]], name = "colocalization")
    if (length(cl) == 0) return(cl)
    cl <- cl %>% left_join(name_conversion[[i]], by = "component") %>%
      select(-component, -sample) %>%
      rename(component = new)
    return(cl)
  }))

  # Convert to CellGraphAssay and add cellgraphs
  pxcellassay <- as.CellGraphAssay(
    x = new_assay,
    cellgraphs = cellgraphs_new,
    polarization = polarization,
    colocalization = colocalization,
    arrow_dir = new_dir,
    arrow_data = arrow_data
  )

  return(pxcellassay)
}

#' Validate polarization \code{tbl_df}
#'
#' @param polarization A \code{tbl_df} with polarization scores
#' @param cell_ids A character vector with cell IDs
#' @param markers A character vector with marker names
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @noRd
.validate_polarization <- function (
  polarization,
  cell_ids,
  markers,
  verbose = FALSE
) {
  # Set polarization to empty tibble if NULL
  polarization <- polarization %||% tibble()
  stopifnot("'polarization' must be a non-empty 'tbl_df' object" = inherits(polarization, "tbl_df"))

  # Validate polarization and colocalization
  if (length(polarization) > 0) {
    # Check column names
    stopifnot("'polarization' names are invalid" =
                all(names(polarization) ==
                      c("morans_i","morans_p_value","morans_p_adjusted","morans_z","marker","component")))
    # Check component names
    cells_in_polarization <- cell_ids %in% (polarization$component %>% unique())
    if (!all(cells_in_polarization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Cells {paste(which(!cells_in_polarization), collapse=', ')} in 'counts' are missing from 'polarization' table")
      return(polarization)
    }
    # Check marker names
    markers_in_polarization <- markers %in% (polarization$marker %>% unique())
    if (!all(markers_in_polarization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Markers {paste(which(!markers_in_polarization), collapse=', ')} in 'counts' are missing from 'polarization' table")
    }
  }
  return(polarization)
}

#' Validate colocalization \code{tbl_df}
#'
#' @param colocalization A \code{tbl_df} with colocalization scores
#' @param cell_ids A character vector with cell IDs
#' @param markers A character vector with marker names
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @noRd
.validate_colocalization <- function (
  colocalization,
  cell_ids,
  markers,
  verbose = FALSE
) {

  # Set colocalization to empty tibble if NULL
  colocalization <- colocalization %||% tibble()
  stopifnot("'colocalization' must be a non-empty 'tbl_df' object" = inherits(colocalization, "tbl_df"))

  if (length(colocalization) > 0) {
    # Check column names
    stopifnot("'colocalization' names are invalid" =
                all(names(colocalization) ==
                      c("marker_1","marker_2","pearson","pearson_mean","pearson_stdev","pearson_z","pearson_p_value",
                        "pearson_p_value_adjusted","jaccard","jaccard_mean","jaccard_stdev","jaccard_z","jaccard_p_value",
                        "jaccard_p_value_adjusted","component")))
    # Check component names
    cells_in_colocalization <- cell_ids %in% (colocalization$component %>% unique())
    if (!all(cells_in_colocalization)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Cells {paste(which(!cells_in_colocalization), collapse=', ')} in 'counts' are missing from 'colocalization' table")
      return(colocalization)
    }
    # Check marker names
    all_markers <- c(colocalization$marker_1 %>% unique(), colocalization$marker_2 %>% unique()) %>% unique()
    if (!all(markers %in% all_markers)) {
      if (verbose && check_global_verbosity())
        cli_alert_warning("Markers {paste(which(!markers %in% all_markers), collapse=', ')} in 'counts' are missing from 'colocalization' table")
    }
  }
  return(colocalization)
}
