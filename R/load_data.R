#' Read a count matrix from a pxl file
#'
#' @param filename Path to a .pxl file
#' @param return_list If TRUE, returns a list with the expression matrix and
#' the output from \code{h5read}
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @return A count matrix or a list if \code{return_list = TRUE}
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' counts <- ReadMPX_counts(pxl_file)
#' counts[1:5, 1:5]
#'
#' @export
#'
ReadMPX_counts <- function(
  filename,
  return_list = FALSE,
  verbose = TRUE
) {
  stopifnot(
    "filename must be a character of length 1" =
      is.character(filename) &&
        (length(filename) == 1)
  )
  if (!file.exists(filename)) abort(glue("{filename} doesn't exist"))

  # Reads anndata from a .pxl file by unzipping it into a temporary folder and reading the anndata object within.
  if (verbose && check_global_verbosity()) {
    cli_alert_info("Loading count data from {filename}")
  }

  # Unzip pxl file
  if (endsWith(filename, ".pxl")) {
    res <- tryCatch(unzip(filename, exdir = fs::path_temp()),
      error = function(e) e,
      warning = function(w) w
    )
    if (inherits(x = res, what = "simpleWarning")) {
      abort("Failed to unzip 'adata.h5ad' data")
    }
  } else {
    abort(glue("Invalid file format .{.file_ext(filename)}. Expected a .pxl file."))
  }

  # Read temporary file
  adata_file <- file.path(fs::path_temp(), "adata.h5ad")
  tmp_file <- fs::file_temp(ext = "h5ad")
  fs::file_move(adata_file, tmp_file)
  hd5_object <- hdf5r::H5File$new(tmp_file, "r")

  # Extract contents
  X <- hd5_object[["X"]]$read()
  colnames(X) <- hd5_object[["obs"]][["component"]]$read()
  markers <- try(
    {
      hd5_object[["var"]][["marker"]]$read()
    },
    silent = TRUE
  )
  if (inherits(markers, "try-error")) {
    markers <- hd5_object[["var"]][["_index"]]$read()
  }
  rownames(X) <- markers

  X <- X[seq_len(nrow(X)), seq_len(ncol(X))]

  hd5_object$close_all()

  if (return_list) {
    return(list(X = X, tmp_file = tmp_file))
  } else {
    return(X)
  }
}

#' Load data from PXL file into a \code{Seurat} object
#'
#' This wrapper function can be used to load data from a PXL file, and returns
#' a \code{Seurat} object.
#'
#' By default, the MPX count matrix is returned in a \code{CellGraphAssay} object.
#' Graphs are not loaded directly unless \code{load_cell_graphs = TRUE}. Graphs
#' can also be loaded at a later stage with \code{\link{LoadCellGraphs}}.
#'
#' When setting the global option \code{Seurat.object.assay.version} to \code{"v5"},
#' the function will return a \code{CellGraphAssay5} object instead.
#'
#' @param assay Assay name
#' @param return_cellgraphassay Should data be loaded as a \code{CellGraphAssay} object?
#' @param load_cell_graphs Should the cellgraphs be loaded into the \code{CellGraphAssay} object?
#' @param load_polarity_scores,load_colocalization_scores Logical specifying if the polarity
#' and colocalization scores should be loaded. These parameters only have an effect if
#' \code{return_cellgraphassay = TRUE}.
#' @param add_additional_assays If other matrix representations are stored in the PXL file,
#' for instance CLR-normalized counts or denoised, set this parameter to \code{TRUE} to load these
#' in separate Assays.
#' @param edgelist_outdir A directory where the edgelist should be stored
#' @param overwrite Should \code{edgelist_outdir} be overwritten?
#' @param ... Additional parameters passed to \code{\link[SeuratObject]{CreateSeuratObject}}
#' @inheritParams ReadMPX_counts
#'
#' @import rlang
#'
#' @family data-loaders
#'
#' @return An object of class \code{Seurat}
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data as a Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' seur_obj <- ReadMPX_Seurat(pxl_file)
#' seur_obj
#'
#' @export
#'
ReadMPX_Seurat <- function(
  filename,
  assay = "mpxCells",
  return_cellgraphassay = TRUE,
  load_cell_graphs = FALSE,
  load_polarity_scores = TRUE,
  load_colocalization_scores = TRUE,
  add_additional_assays = FALSE,
  edgelist_outdir = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  stopifnot(
    "assay must be a character of length 1" =
      is.character(assay) &&
        (length(assay) == 1)
  )

  # Load count matrix
  data <- ReadMPX_counts(filename = filename, return_list = TRUE, verbose = FALSE)
  X <- data$X
  hd5_object <- hdf5r::H5File$new(data$tmp_file, "r")

  # Load edgelist
  empty_graphs <- rep(list(NULL), ncol(X)) %>% set_names(nm = colnames(X))

  if (return_cellgraphassay) {
    # Create CellGraphAssay(5)
    cg_assay_create_func <- switch(getOption("Seurat.object.assay.version", "v3"),
      "v3" = CreateCellGraphAssay,
      "v5" = CreateCellGraphAssay5
    )
    cg_assay <- cg_assay_create_func(
      counts = X,
      cellgraphs = empty_graphs,
      fs_map = tibble(
        id_map = list(tibble(
          current_id = colnames(X),
          original_id = colnames(X)
        )),
        sample = 1L,
        pxl_file = ifelse(.is_absolute_path(filename), normalizePath(filename), filename)
      )
    )
    Key(cg_assay) <- paste0(assay, "_")

    if (load_cell_graphs) {
      if (verbose && check_global_verbosity()) {
        cli_alert_warning(
          glue(col_br_red(
            "Loading cell graphs into the Seurat object may take time and uses a lot of memory.",
            " Consider using LoadCellGraphs to load selected cell graphs at a later stage."
          ))
        )
      }
      cg_assay <- LoadCellGraphs(cg_assay, verbose = verbose)
    }
    # Load polarity scores
    if (load_polarity_scores) {
      polarization <- ReadMPX_item(filename = filename, items = "polarization", verbose = FALSE)
      if (inherits(polarization$component, "integer")) {
        polarization$component <- as.character(polarization$component)
      }
      cg_assay@polarization <- polarization
    }
    # Load colocalization scores
    if (load_colocalization_scores) {
      colocalization <- ReadMPX_item(filename = filename, items = "colocalization", verbose = FALSE)
      if (inherits(colocalization$component, "integer")) {
        colocalization$component <- as.character(colocalization$component)
      }
      # TODO: Remove this once the colocalization tables have been updated
      if (all(c("marker1", "marker2") %in% names(colocalization))) {
        colocalization <- colocalization %>%
          rename(marker_1 = !!sym("marker1"), marker_2 = !!sym("marker2"))
      }
      cg_assay@colocalization <- colocalization
    }
  } else {
    cg_assay <- switch(getOption("Seurat.object.assay.version", "v3"),
      "v3" = CreateAssayObject(counts = X),
      "v5" = CreateAssay5Object(counts = X)
    )
    Key(cg_assay) <- paste0(assay, "_")
  }

  # Create Seurat object
  seur_obj <- CreateSeuratObject(counts = cg_assay, assay = assay, ...)
  if (verbose && check_global_verbosity()) {
    cli_alert_success(glue(
      "Created a 'Seurat' object with {col_br_blue(ncol(X))} cells ",
      "and {col_br_blue(nrow(X))} targeted surface proteins"
    ))
  }

  # Extract meta data
  meta_data <-
    names(hd5_object[["obs"]]) %>%
    lapply(function(nm) {
      if (length(names(hd5_object[["obs"]][[nm]])) == 2) {
        indices <- hd5_object[["obs"]][[nm]][["codes"]]$read() + 1
        categories <- hd5_object[["obs"]][[nm]][["categories"]]$read()
        col <- tibble(!!sym(nm) := factor(categories[indices], categories))
      } else {
        col <- tibble(!!sym(nm) := hd5_object[["obs"]][[nm]]$read())
      }
      return(col)
    }) %>%
    do.call(bind_cols, .) %>%
    as.data.frame() %>%
    column_to_rownames("component")

  if (!all(rownames(meta_data) == colnames(seur_obj))) {
    abort(glue(
      "The cell names in the Seurat object and the metadata do not match. ",
      "\nFailed to create Seurat object."
    ))
  }
  seur_obj@meta.data <- meta_data

  # Extract feature meta data (var)
  feature_meta_data <-
    names(hd5_object[["var"]]) %>%
    lapply(function(nm) {
      # If reading function is missing, don't read the column
      # This is likely due to the column being empty
      col <- try({tibble(!!sym(nm) := hd5_object[["var"]][[nm]]$read())}, silent = TRUE)
      if (inherits(col, "try-error") {
        warn(glue("Column '{nm}' in var is empty. Skipping."))
        return(NULL)
      }
      return(col)
    }) %>%
    do.call(bind_cols, .) %>%
    as.data.frame()

  # Close hdf5 file
  hd5_object$close()

  if (getOption("Seurat.object.assay.version", "v3") == "v3") {
    rownames(feature_meta_data) <- feature_meta_data$marker
    seur_obj[[assay]]@meta.features <- feature_meta_data
  } else {
    seur_obj[[assay]]@meta.data <- feature_meta_data
  }

  # Set variable features to all features
  VariableFeatures(seur_obj) <- rownames(seur_obj)

  # Remove redundant meta data columns
  seur_obj@meta.data <- seur_obj@meta.data %>% select(-contains(c("nFeature", "nCount", "_index")))

  return(seur_obj)
}

#' Read an MPX data item
#'
#' \code{ReadMPX_item} reads any number of items from a .pxl file. If multiple
#' items are specified, the output is a list of  \code{tbl_df} objects.
#' Otherwise the output is a single \code{tbl_df}. \code{ReadMPX_polarization},
#' \code{ReadMPX_colocalization} and \code{ReadMPX_edgelist} are wrappers for
#' \code{ReadMPX_item} to read specific items.
#'
#' @param filename Path to a .pxl file
#' @param items One or several of "colocalization", "polarization", "edgelist"
#' @param verbose Print messages
#'
#' @import rlang
#'
#' @family data-loaders
#'
#' @rdname ReadMPX_item
#'
#' @return Either an object of class \code{tbl_df} or a list of \code{tbl_df}
#' objects if multiple \code{items} are selected
#'
#' @examples
#'
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' polarization <- ReadMPX_item(pxl_file, items = "polarization")
#' polarization
#'
#' # Alternative 2
#' polarization <- ReadMPX_polarization(pxl_file)
#'
#' @export
#'
ReadMPX_item <- function(
  filename,
  items = c("colocalization", "polarization", "edgelist"),
  verbose = TRUE
) {
  # Check input parameters
  stopifnot(
    "filename must be a character of length 1" =
      is.character(filename) &&
        (length(filename) == 1),
    "Invalid items" =
      all(items %in% c("colocalization", "polarization", "edgelist")),
    "Expected a .pxl file" =
      .file_ext(filename) == "pxl"
  )
  if (!file.exists(filename)) abort(glue("{filename} doesn't exist"))

  if (verbose && check_global_verbosity()) {
    cli_alert_info("Loading item(s) from: {filename}")
  }

  suppressWarnings(
    read_items <-
      items %>%
      set_names(., .) %>%
      lapply(function(item) {
        if (verbose && check_global_verbosity()) {
          cli_alert("  Loading {col_br_magenta(item)} data")
        }

        # generate a file name for a new temporary directory.
        # When running R CMD check on windows latest release, we
        # do not have permissions to modify files in TEMPDIR.
        # For this reason, we create a new directory with a unique
        # name instead.
        exdir_temp <- fs::file_temp()

        item_name <- paste0(item, ".parquet")

        # Unzip item to temporary directory
        unzipped_filename <- unzip(filename, item_name, exdir = exdir_temp)

        if (length(unzipped_filename) == 0) {
          abort(glue("Failed to extract {item_name} from {filename}"))
        }

        # Read contents of parquet file
        outdata <- read_parquet(unzipped_filename)

        # Try to delete temporary directory or throw a warning if it fails.
        check <- tryCatch(
          {
            fs::dir_delete(exdir_temp)
          },
          error = function(e) {
            e
          },
          warning = function(w) {
            w
          }
        )
        if (inherits(check, what = "error")) {
          cli_alert_warning("Failed to remove temporary directory {exdir_temp}")
        }

        return(outdata)
      })
  )

  if (length(read_items) == 1) {
    read_items <- read_items[[1]]
  }

  if (verbose && check_global_verbosity()) {
    if (inherits(read_items, what = "list")) {
      cli_alert_success("Returning a list of tibbles with items: {paste(items, collapse = ', ')}")
    } else {
      cli_alert_success("Returning a '{class(read_items)[1]}' object")
    }
  }
  return(read_items)
}


#' @rdname ReadMPX_item
#'
#' @export
#'
ReadMPX_polarization <- function(
  filename,
  verbose = TRUE
) {
  ReadMPX_item(filename, items = "polarization", verbose = verbose)
}

#' @rdname ReadMPX_item
#'
#' @export
#'
ReadMPX_colocalization <- function(
  filename,
  verbose = TRUE
) {
  ReadMPX_item(filename, items = "colocalization", verbose = verbose)
}

#' @rdname ReadMPX_item
#'
#' @export
#'
ReadMPX_edgelist <- function(
  filename,
  verbose = TRUE
) {
  ReadMPX_item(filename, items = "edgelist", verbose = verbose)
}

#' Read metadata from a PXL file
#'
#' @param filename Path to a PXL file
#'
#' @rdname ReadMPX_metadata
#'
#' @examples
#' library(pixelatorR)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' meta_data <- ReadMPX_metadata(pxl_file)
#'
#' # Check pixelator version and sample ID
#' meta_data
#'
#' # Check parameter settings for the pixelator run
#' meta_data$analysis$params
#'
#' @export
#'
ReadMPX_metadata <- function(
  filename
) {
  if (!fs::file_exists(filename)) {
    abort(glue("File {col_blue(filename)} doesn't exist"))
  }

  if (fs::path_ext(filename) != "pxl") {
    abort(glue("File {col_blue(filename)} is not a PXL file"))
  }

  temp_dir <- fs::path_temp()
  temp_file <- fs::path(temp_dir, "metadata.json")
  zip::unzip(
    zipfile = filename,
    files = "metadata.json",
    exdir = temp_dir
  )

  # Read meta data from JSON file
  meta_data <- jsonlite::read_json(temp_file, simplifyVector = TRUE) %>%
    as_tibble() %>%
    select(-contains("file_format_version"))
  if (nrow(meta_data) == 0) {
    abort("No metadata found in the PXL file")
  }
  analysis <- meta_data$analysis[[1]]
  if (any(c("polarization", "colocalization") %in% names(analysis))) {
    analysis <- unlist(analysis %>% unname())
    names(analysis) <- names(analysis) %>%
      stringr::str_replace("\\.", "_")
  } else {
    analysis <- unlist(analysis)
  }
  meta_data$analysis <- list(params = analysis)
  class(meta_data) <- c("pixelator_metadata", class(meta_data))

  return(meta_data)
}


#' Print method for pixelator_metadata
#'
#' @param x An object of class \code{pixelator_metadata}
#' @param detailed Logical. If \code{TRUE} print all metadata, if
#' \code{FALSE} print only the pixelator version and sample ID.
#' @param ... Additional arguments passed to \code{\link{print}}
#'
#' @rdname print.pixelator_metadata
#' @docType methods
#'
#' @examples
#' library(pixelatorR)
#' library(dplyr)
#'
#' # Load example data
#' pxl_file <- system.file("extdata/five_cells",
#'   "five_cells.pxl",
#'   package = "pixelatorR"
#' )
#' meta_data <- ReadMPX_metadata(pxl_file)
#'
#' # Check pixelator version and sample ID
#' meta_data
#'
#' # Multiple files with less detail
#' pxl_files <- c(pxl_file, pxl_file)
#' meta_data_merged <- lapply(pxl_files, ReadMPX_metadata) %>%
#'   bind_rows()
#' meta_data_merged %>% print(detailed = FALSE)
#'
#' @export
#'
print.pixelator_metadata <- function(
  x,
  detailed = TRUE,
  ...
) {
  sampleIDs <- x$sample
  pixelator_version <- x$version

  if ("params" %in% names(x$analysis)) {
    analysis_params <- lapply(x$analysis, function(x) {
      x %>% unlist()
    })
  } else {
    analysis_params <- NULL
  }

  for (i in seq_len(nrow(x))) {
    glue(
      "Sample {i} name: \t\t{col_blue(sampleIDs[i])}\n",
      "Pixelator version: \t{col_br_magenta(pixelator_version[i])}\n",
      .trim = FALSE
    ) %>% print()

    if (!is.null(analysis_params) && detailed) {
      cli_h2("Analysis parameters")
      analysis_params_cur <- tibble(
        parameter = names(analysis_params[[i]]),
        value = analysis_params[[i]], row.names = NULL
      )
      analysis_params_cur %>% print()
    }

    if (length(analysis_params) > 1 && i < length(analysis_params)) {
      if (detailed) {
        cat("\n")
      }
      cli::cat_rule()
    }
  }
}
