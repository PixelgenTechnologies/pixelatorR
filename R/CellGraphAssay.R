#' @include generics.R
#' @importFrom methods setClass setClassUnion setMethod slot slot<- new as slotNames
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
    stopifnot(
      "column 'component' is missing from cellgraphs" =
        "component" %in% names(fsd)
    )
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
# Get/set methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
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
CellGraphs.CellGraphAssay <- function (
  object,
  ...
) {
  return(slot(object, name = "cellgraphs"))
}


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
      inherits(new.names, what = "character") &&
      (length(new.names) == ncol(object))
  )

  # save original cell IDs
  orig.names <- colnames(object)

  # Fetch unique sample IDs from orig.names and new.names
  new_sample_id_table <- do.call(rbind, strsplit(new.names, "_"))
  new_sample_id <- new_sample_id_table[, 1] %>% unique()
  old_sample_id_table <- do.call(rbind, strsplit(orig.names, "_"))
  if (ncol(old_sample_id_table) == 1) {
    old_sample_id <- "S1"
  } else {
    old_sample_id <- old_sample_id_table[, 1] %>% unique()
  }

  # Validate names
  names_checked <- sapply(new.names, function(s) {
    stringr::str_like(s, pattern = "^[a-zA-Z][a-zA-Z0-9]*\\_RCVCMP\\d{7}$")
  })
  if (!any(names_checked)) {
    abort(glue("Failed to merge CellGraphAssays.\n\n",
               "Make sure to follow these steps:\n\n",
               "1. CellGraphAssay column names should have the following format:\n",
               "   ^[a-zA-Z][a-zA-Z0-9]*_RCVCMP\\d{7}$\n",
               "   where the first part is a sample ID and the second part is the PXL ID, \n",
               "   separated by an underscore. For example, {col_green('Sample1_RCVCMP0000000')}.\n\n",
               "2. When merging Seurat objects, make sure to set {col_green('add.cell.ids')}.\n\n",
               "3. Cannot merge Seurat objects with CellGraphAssays twice due to naming conflicts. \n",
               "   Instead, merge all data sets once:\n",
               "   {col_green('se_merged <- merge(se1, list(se2, se3, se4, ...))')}\n",
               "   Attempting to merge the merged object again will fail:\n",
               "   {col_red('se_double_merged <- merge(se_merged, se_merged)')}"))
  }

  names(slot(object = object, name = "cellgraphs")) <- new.names

  names(x = new.names) <- NULL
  for (data.slot in object[]) {
    old.data <- LayerData(object = object, layer = data.slot)
    if (ncol(x = old.data) <= 1) {
      next
    }
    colnames(x = slot(object = object, name = data.slot)) <- new.names
  }

  # Get arrow dir
  arrow_dir <- ArrowDir(object)
  if (!is.null(arrow_dir)) {
    arrow_dirs <- list.files(arrow_dir, full.names = TRUE)

    # Create new directory
    session_tmpdir_random <-
      file.path(
        getOption("pixelatorR.arrow_outdir"),
        glue("{.generate_random_string()}-{format(Sys.time(), '%Y-%m-%d-%H%M%S')}"))

    # Trigger garbage cleaning if the edgelist directories exceed the
    # maximum allowed size
    .run_clean()

    if (length(arrow_dirs) == 1) {
      # Handle renaming if 1 hive-style directory is present

      # Check sample ID
      if (length(new_sample_id) > 1) {
        abort(glue("Found multiple sample IDs in 'new.names' but only 1 edgelist in the arrow directory."))
      }

      # Copy directory
      if (dir.exists(arrow_dirs)) {
        dir.create(session_tmpdir_random)
        file.copy(from = arrow_dirs, to = session_tmpdir_random, recursive = TRUE)
        # Rename hive-style directory
        hive_style_dir_sample1 <- list.files(session_tmpdir_random, full.names = TRUE)
        file.rename(from = hive_style_dir_sample1, file.path(session_tmpdir_random, paste0("sample=", new_sample_id)))
      } else {
        abort(glue("Directory '{arrow_dirs}' is missing Cannot rename cell IDs in edgelists."))
      }
    } else {
      # Handle renaming if more than 1 hive-style directories are present

      if (!length(arrow_dirs) == length(new_sample_id)) {
        abort(glue("Found {length(arrow_dirs)} samples in arrow directory, but {length(new_sample_id)} samples in 'new.names'"))
      }

      # Copy directories to new folder
      if (all(dir.exists(arrow_dirs))) {
        dir.create(session_tmpdir_random)
        for (i in seq_along(arrow_dirs)) {
          file.copy(from = arrow_dirs[i], to = session_tmpdir_random, recursive = TRUE)
        }
        hive_style_dir_samples <- list.files(session_tmpdir_random, full.names = TRUE)
        hive_style_dir_sample_IDs <- basename(hive_style_dir_samples) %>% gsub(pattern = "sample=", replacement = "", x = .)
        hive_style_dir_samples <- setNames(hive_style_dir_samples, nm = hive_style_dir_sample_IDs)

        # Reorder hive_style_dir_samples
        hive_style_dir_samples <- hive_style_dir_samples[old_sample_id]

        # Rename hive-style directory
        for (i in seq_along(arrow_dirs)) {
          file.rename(from = hive_style_dir_samples[i], file.path(session_tmpdir_random, paste0("sample=", new_sample_id[i])))
        }
      } else {
        abort(glue("The following directories are missing:\n {paste0(arrow_dirs, collapse='\n')} ",
                   "\nCannot rename cell IDs in edgelists."))
      }
    }

    # Log command
    command <- sys.calls()[[1]][1] %>% as.character()
    options(pixelatorR.edgelist_copies = bind_rows(
      getOption("pixelatorR.edgelist_copies"),
      tibble(command = command,
             edgelist_dir = normalizePath(session_tmpdir_random),
             timestamp = Sys.time())
    ))

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
  arrow_dir = NULL,
  arrow_data = NULL,
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

  # Validate polarization and colocalization
  polarization <- .validate_polarization(polarization, cell_ids = colnames(x), markers = rownames(x))
  colocalization <- .validate_colocalization(colocalization, cell_ids = colnames(x), markers = rownames(x))

  # Abort if cellgraphs is empty and neither arrow_dir or arrow_data is provided
  loaded_graphs <- sum(sapply(cellgraphs, is.null))
  if (loaded_graphs == ncol(x)) {
    stopifnot(
      "One of 'arrow_dir' or 'arrow_data' must be provided if 'cellgraphs is empty'" =
        (!is.null(arrow_dir)) ||
        (!is.null(arrow_data))
    )
  }

  # Handle arrow_dir
  arrow_dir <- arrow_dir %||% NA_character_
  if (!is.na(arrow_dir)) {
    stopifnot(
      "'arrow_dir' must be a non-empty character" =
        is.character(arrow_dir) &&
        (length(arrow_dir) >= 1)
    )
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
    #stopifnot("One or several components are missing from 'arrow_data'" = all(colnames(x) %in% (arrow_data %>% pull(component, as_vector = TRUE))))
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
#' # Set arrow data output directory to temp for tests
#' options(pixelatorR.arrow_outdir = tempdir())
#'
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
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
#' seur_subset <- subset(seur, cells = colnames(seur)[1:3])
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

  stopifnot(
    "All 'cells' must be present in x" =
      all(cells %in% colnames(x))
  )

  # Get cellgraphs
  cellgraphs <- x@cellgraphs

  # subset elements in the standard assay
  standardassay <- as(object = x, Class = "Assay")
  standardassay <- subset(x = standardassay, features = features, cells = cells)

  # Filter cellgraphs
  cellgraphs_filtered <- cellgraphs[colnames(standardassay)]

  # Fetch arrow_dir
  arrow_dir <- slot(x, name = "arrow_dir")

  # Get sample ids
  sample_id_table <- do.call(rbind, strsplit(colnames(x), "_"))
  rownames(sample_id_table) <- colnames(x)
  if (ncol(sample_id_table) == 1) {
    sample_id <- NULL
  } else {
    sample_id <- sample_id_table[, 1] %>% unique()
  }

  # Handle arrow data set if available
  if (!is.na(arrow_dir)) {
    x <- RestoreArrowConnection(x, verbose = FALSE)

    # Filter edgelists
    if (!is.null(cells)) {

      # Only filter by cells
      if (length(cells) < ncol(x)) {

        # Trigger garbage cleaning if the edgelist directories exceed the
        # maximum allowed size
        .run_clean()

        # Create a temporary directory with a unique name
        session_tmpdir_random <-
          file.path(
            getOption("pixelatorR.arrow_outdir"),
            glue("{.generate_random_string()}-{format(Sys.time(), '%Y-%m-%d-%H%M%S')}"))
        dir.create(session_tmpdir_random)

        # Handle samples
        if (ncol(sample_id_table) > 1) {
          components_keep_list <- split(sample_id_table[cells, 2], sample_id_table[cells, 1])
          for (s in sample_id) {
            components_keep <- components_keep_list[[s]]
            slot(x, name = "arrow_data") %>%
              filter(sample == s) %>%
              filter(component %in% components_keep) %>%
              group_by(sample) %>%
              write_dataset(path = session_tmpdir_random)
          }
        } else {
          # Filter edgelist and export it
          slot(x, name = "arrow_data") %>%
            filter(component %in% cells) %>%
            group_by(sample) %>%
            write_dataset(session_tmpdir_random)
        }

        # Rename parquet files for consistensy with other functions
        files <- list.files(session_tmpdir_random, pattern = "parquet", recursive = TRUE, full.names = TRUE)
        for (f in files) {
          if (basename(f) != "edgelist.parquet") {
            file.rename(from = f, file.path(dirname(f), "edgelist.parquet"))
          }
        }

        # Log command
        command <- sys.calls()[[1]][1] %>% as.character()
        options(pixelatorR.edgelist_copies = bind_rows(
          getOption("pixelatorR.edgelist_copies"),
          tibble(command = command,
                 edgelist_dir = normalizePath(session_tmpdir_random),
                 timestamp = Sys.time())
        ))

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


#' @param add.cell.ids A character vector with sample names
#'
#' @importFrom SeuratObject Key Key<- RenameCells
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
    add.cell.ids = NULL,
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

  # Define add.cell.ids
  # add.cell.ids <- add.cell.ids %||% paste0("Sample", seq_along(objects))
  if (!is.null(add.cell.ids)) {
    stopifnot(
      "Length of 'add.cell.ids' must match the number of objects to merge" =
        length(add.cell.ids) == length(objects)
    )
  }

  # Check duplicate cell names
  cell.names <- unlist(lapply(objects, colnames))
  unique_names <- table(cell.names)
  names_are_duplicated <- any(unique_names > 1)
  sample_id_old_table <- do.call(rbind, strsplit(cell.names, "_"))
  if (names_are_duplicated && is.null(add.cell.ids)) {
    if (ncol(sample_id_old_table) == 1) {
      add.cell.ids <- paste0("Sample", seq_along(objects))
    } else if (ncol(sample_id_old_table) == 2) {
      abort("Found non-unique IDs across samples. A 'add.cell.ids' must be specified.")
    }
  }

  # Fetch sample IDs from column names
  if (ncol(sample_id_old_table) == 1) {
    objects <-
      lapply(seq_along(objects), function(i) {
        if (is.null(add.cell.ids)) {
          return(objects[[i]])
        } else {
          new_names_modified <- paste0(add.cell.ids[i], "_", Cells(x = objects[[i]]))
          return(RenameCells(object = objects[[i]], new.names = new_names_modified))
        }
      })
  } else {
    objects <-
      lapply(seq_along(objects), function(i) {
        if (is.null(add.cell.ids)) {
          return(objects[[i]])
        } else {
          cli_alert_warning("Found multiple samples in objects. 'add.cell.ids' will be added as a prefix to old IDs.")
          new_names_modified <- paste0(add.cell.ids[i], Cells(x = objects[[i]]))
          return(RenameCells(object = objects[[i]], new.names = new_names_modified))
        }
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

    # Trigger garbage cleaning if the edgelist directories exceed the
    # maximum allowed size
    .run_clean()

    if (!any(dir.exists(all_arrow_dirs)))
      abort(glue("Paths to edgelist parquet file directories {paste(all_arrow_dirs, collapse = ',')} doesn't exist"))

    # Create a new temporary directory with a time stamp
    new_dir <-
      file.path(
        getOption("pixelatorR.arrow_outdir"),
        glue("{.generate_random_string()}-{format(Sys.time(), '%Y-%m-%d-%H%M%S')}"))
    dir.create(path = new_dir, showWarnings = FALSE)

    # Move hive-style old sample diretories to new directory
    for (i in seq_along(all_arrow_dirs)) {
      hive_style_dirs <- list.files(all_arrow_dirs[i], full.names = TRUE)
      for (ii in seq_along(hive_style_dirs)) {
        file.copy(from = hive_style_dirs[ii], to = new_dir, recursive = TRUE)
      }
    }

    # Log command
    command <- sys.calls()[[1]][1] %>% as.character()
    options(pixelatorR.edgelist_copies = bind_rows(
      getOption("pixelatorR.edgelist_copies"),
      tibble(command = command,
             edgelist_dir = normalizePath(new_dir),
             timestamp = Sys.time())
    ))

    # Open arrow data
    arrow_data <- open_dataset(new_dir)
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
