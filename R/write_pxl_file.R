# Declarations used in package check
globalVariables(
  names = c('original_id', 'current_id'),
  package = 'pixelatorR',
  add = TRUE
)

#' Export Seurat object data to a .pxl file
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The function only exports the essential data required to create a functional
#' .pxl file. See details below for a description of what data is exported.
#'
#' @section Exported data:
#'
#'  - Count data and metadata : The raw count matrix and component metadata from the Seurat object
#'  are exported into a .h5ad file that can be read with the anndata Python library.
#'  - Polarization scores : The polarization scores are exported to a .parquet file.
#'  - Colocalization scores : The colocalization scores are exported to a .parquet file.
#'  - Edgelist data : The edgelist data is first collected from the original .pxl file(s),
#'  filtered to include the components currently available in the \code{Seurat} object,
#'  and then the component IDs are updated. The resulting merged edgelist data is exported
#'  to a .parquet file.
#'  - Sample meta data : The sample meta data is extracted from the original .pxl file(s),
#'  then merged and exported to a .json file.
#'
#' The structure of the .pxl file is detailed below:
#'
#' |-- adata.h5ad\cr
#' |-- polarization.parquet\cr
#' |-- colocalization.parquet\cr
#' |-- metadata.json\cr
#' |-- edgelist.parquet
#'
#' The merged files are converted into a zip archive and saved to the target .pxl file.
#'
#' @param object A \code{Seurat} object with a \code{CellGraphAssay5}
#' assay object created with pixelatorR.
#' @param file A character string specifying the path to the .pxl
#' file to be created.
#' @param assay A character string specifying the name of the \code{CellGraphAssay5}.
#' If set to \code{NULL} the default assay will be used.
#' @param export_layouts A logical value specifying whether to export the layouts from
#' the \code{Seurat} object. If set to \code{TRUE}, each component must have a layout.
#' See \code{\link{LoadCellGraphs}} and \code{\link{ComputeLayout}} for details about
#' how to load graphs and compute layouts.
#' @param overwrite A logical value specifying whether to overwrite the \code{file}
#' if it already exists.
#'
#' @return Nothing. The function writes the .pxl file to the specified location.
#'
#' @examples
#' # Use Assay5 as the default assay version
#' options(Seurat.object.assay.version = "v5")
#'
#' # Create Seurat object
#' pxl_file <- system.file("extdata/five_cells",
#'                         "five_cells.pxl",
#'                         package = "pixelatorR")
#' se <- ReadMPX_Seurat(pxl_file)
#'
#' se_merged <- merge(se, list(se, se, se))
#' pxl_file <- fs::path_temp("small.pxl")
#'
#' # Export data to a new .pxl file
#' WriteMPX_pxl_file(se_merged, pxl_file)
#'
#' # Read the new .pxl file
#' se_merged <- ReadMPX_Seurat(pxl_file)
#'
#' # Reset global option
#' options(Seurat.object.assay.version = "v3")
#'
#' @export
#'
WriteMPX_pxl_file <- function (
  object,
  file,
  assay = NULL,
  export_layouts = FALSE,
  overwrite = FALSE
) {

  stopifnot(
    "'object' must be a Seurat object" =
      inherits(object, what = "Seurat"),
    "'file' must be a valid path ending with .pxl" =
      is.character(file) && (fs::path_ext(file) == "pxl"),
    "'overwrite' must be either TRUE or FALSE" =
      is.logical(overwrite)
  )

  # Check if file exists
  if (fs::file_exists(file) & !overwrite) {
    abort(glue("{col_br_blue(file)} already exists. Please select a different ",
               "file name or set overwrite = TRUE if you are certain that the ",
               "existing file should be replaced."))
  }

  assay <- assay %||% DefaultAssay(object)
  cg_assay <- object[[assay]]
  stopifnot(
    "'assay' must be a 'CellGraphAssay' or a 'CellGraphAssay5' object" =
      is(cg_assay, "MPXAssay")
  )

  # fetch and validate fs_map
  fs_map <- cg_assay@fs_map
  .validate_fs_map(fs_map)

  # Delete the file only if it's not in use
  # by the current object
  if (fs::file_exists(file)) {
    if (file %in% fs_map$pxl_file) {
      abort(glue(
        "The selected file name '{col_br_blue(file)}' is currently in use ",
        "by the input object. You need to select a different file name."
      ))
    }
    fs::file_delete(file)
  }

  # Create a temporary directory
  pxl_folder <- file.path(fs::path_temp(), .generate_random_string())
  fs::dir_create(path = pxl_folder)

  cli_rule(left = "Creating .pxl file")

  # Export spatial metrics
  .write_spatial_metrics_to_parquet(object = object, pxl_folder = pxl_folder)

  # Update and write edgelist data
  .merge_and_write_edgelists_to_parquet(fs_map = fs_map, pxl_folder = pxl_folder)

  # Extract and merge meta data
  .write_json_meta_data(fs_map = fs_map, pxl_folder = pxl_folder)

  # Extract Seurat data to adata.h5ad
  .write_to_h5ad_file(object = object,
                      cg_assay = cg_assay,
                      pxl_folder = pxl_folder)

  # Export layout data
  if (export_layouts) {
    cg_list <- CellGraphs(object)
    graphs_check <- sapply(cg_list, function(cg) !is.null(cg))
    if (!all(graphs_check)) {
      abort(glue("CellGraphs must be available for all {length(cg_list)} components. ",
                 "Found CellGraphs for {sum(graphs_check)} components. \n",
                 "You can either run LoadCellGraphs and ComputeLayout to ",
                 "obtain layouts or set export_layouts = FALSE if you don't want to export the layouts."))
    }
    layouts_check <- sapply(cg_list, function(cg) !is.null(cg@layout))
    if (all(layouts_check)) {
      .merge_layout_with_counts_and_write_to_parquet(cg_list = cg_list, pxl_folder = pxl_folder)
    } else {
      abort(glue("Layouts must be available for all {length(cg_list)} components. ",
                 "Found layouts for {sum(layouts_check)} components. "))
    }
  }

  # Zip the pxl folder and move to the target file
  # We use the simplest compression here since the
  # files are already compressed
  cli_alert_info("Saving .pxl file to {file}")

  zip::zip(
    zipfile = file,
    files = list.files(pxl_folder),
    root = pxl_folder,
    compression_level = 0,
    recurse = TRUE,
    include_directories = FALSE
  )

  cli_alert_success("Finished!")

}

#' Utility function used to fetch the count matrix
#' from a \code{CellGraphAssay} or \code{CellGraphAssay5} object.
#'
#' @param cg_assay A \code{CellGraphAssay} or \code{CellGraphAssay5} object
#'
#' @return A count matrix
#'
#' \code{CellGraphAssay5} objects inherits the \code{Assay5} class which
#' stored the counts matrices in separate layers. This function will
#' combine all layers into a single matrix.
#'
#' @noRd
#'
.fetch_counts <- function (
  cg_assay
) {
  if (inherits(cg_assay, "Assay")) {
    X <- LayerData(cg_assay, layer = "counts")
  } else {
    # Fetch available layers
    all_layers <- Layers(cg_assay)
    # Filter out layers that are not counts (raw count data)
    all_layers <- all_layers[str_detect(all_layers, "^counts\\.[0-9]+$")]
    # Sort layers by their sample index
    all_layers <- all_layers %>% set_names(all_layers)
    all_layers <- all_layers[paste0("counts.", seq_along(all_layers))]
    X <- try({
      lapply(all_layers, function(lr) {
        LayerData(cg_assay, layer = lr)
      }) %>% Reduce(cbind, .)
    }, silent = TRUE)
    if (inherits(X, "try-error"))
      abort(glue("Failed to combine 'Assay5' layers due to invalid dimensions. \n",
                  "Please ensure that all layers have the same number of rows."))
  }
  return(X)
}


#' Utility function used to write Seurat data to a .h5ad file
#'
#' @param object A \code{Seurat} object
#' @param cg_assay A \code{CellGraphAssay} or \code{CellGraphAssay5} object
#' @param pxl_folder A character string specifying the path to a temporary
#' directory where the .h5ad file will be created.
#'
#' @return Nothing
#'
#' @references https://cran.r-project.org/web/packages/hdf5r/vignettes/hdf5r.html
#'
#' @noRd
#'
.write_to_h5ad_file <- function (
  object,
  cg_assay,
  pxl_folder
) {

  # Create a writable .h5ad file
  adata_new <- hdf5r::H5File$new(file.path(pxl_folder, "adata.h5ad"), "w")

  # Fetch raw counts from CellgraphAssay(5) object
  X <- .fetch_counts(cg_assay)
  adata_new$create_dataset("X", robj = as.matrix(X),
                           dims = dim(X),
                           dtype = hdf5r::h5types$double,
                           chunk_dims = NULL)
  hdf5r::h5attr(adata_new[["X"]], which = "encoding-type") <- "array"
  hdf5r::h5attr(adata_new[["X"]], which = "encoding-version") <- "0.2.0"

  # Create layers group (currently empty)
  layers <- adata_new$create_group("layers")
  hdf5r::h5attr(layers, which = "encoding-type") <- "dict"
  hdf5r::h5attr(layers, which = "encoding-version") <- "0.1.0"

  # Create obs group
  obs <- adata_new$create_group("obs")
  obs$create_dataset("component", robj = colnames(object), chunk_dims = NULL)
  hdf5r::h5attr(obs, which = "_index") <- "component"

  # Only export selected meta data columns if they are available
  cols_of_interest <-
    c("vertices", "edges", "molecules", "antibodies", "upia", "upib", "umi", "reads",
      "mean_reads", "median_reads", "mean_upia_degree", "median_upia_degree",
      "mean_umi_per_upia", "median_umi_per_upia", "umi_per_upia", "upia_per_upib",
      "leiden", "tau_type", "tau")
  obs_cols_keep <- intersect(
    cols_of_interest,
    colnames(object[[]])
  )
  hdf5r::h5attr(obs, which = "col-order") <- obs_cols_keep
  if (length(obs_cols_keep) > 0) {
    for (col_name in obs_cols_keep) {
      obs$create_dataset(col_name, robj = object[[]] %>% pull(col_name), chunk_dims = NULL)
    }
  }
  hdf5r::h5attr(obs, which = "encoding-type") <- "dataframe"
  hdf5r::h5attr(obs, which = "encoding-version") <- "0.2.0"

  # Create obsm H5Group (currently empty)
  obsm <- adata_new$create_group("obsm")
  hdf5r::h5attr(obsm, which = "encoding-type") <- "dict"
  hdf5r::h5attr(obsm, which = "encoding-version") <- "0.1.0"

  # Create obsp H5Group (currently empty)
  obsp <- adata_new$create_group("obsp")
  hdf5r::h5attr(obsp, which = "encoding-type") <- "dict"
  hdf5r::h5attr(obsp, which = "encoding-version") <- "0.1.0"

  # Create uns H5Group (currently empty)
  uns <- adata_new$create_group("uns")
  hdf5r::h5attr(uns, which = "encoding-type") <- "dict"
  hdf5r::h5attr(uns, which = "encoding-version") <- "0.1.0"

  # Create var H5Group
  var <- adata_new$create_group("var")
  if (is(cg_assay, "CellGraphAssay5")) {
    feature_meta_data <- cg_assay@meta.data
  } else {
    feature_meta_data <- cg_assay@meta.features
  }

  # Only export selected faeture meta data columns if they are available
  if (length(feature_meta_data) > 0) {
    var$create_dataset("marker", robj = feature_meta_data %>% pull(marker), chunk_dims = NULL)
    var_cols_of_interest <- c("antibody_count", "components", "antibody_pct", "nuclear", "control")
    var_cols_keep <- intersect(
      var_cols_of_interest,
      colnames(feature_meta_data)
    )
    if (length(var_cols_keep) > 0) {
      for (col_name in var_cols_keep) {
        var$create_dataset(col_name, robj = feature_meta_data %>% pull(col_name), chunk_dims = NULL)
      }
    }
  } else {
    var$create_dataset("marker", robj = cg_assay %>% rownames(), chunk_dims = NULL)
  }
  hdf5r::h5attr(var, which = "_index") <- "marker"
  hdf5r::h5attr(var, which = "encoding-type") <- "dataframe"
  hdf5r::h5attr(var, which = "encoding-version") <- "0.2.0"

  # Create varm H5Group (currently empty)
  varm <- adata_new$create_group("varm")
  hdf5r::h5attr(varm, which = "encoding-type") <- "dict"
  hdf5r::h5attr(varm, which = "encoding-version") <- "0.1.0"

  # Create varp H5Group (currently empty)
  varp <- adata_new$create_group("varp")
  hdf5r::h5attr(varp, which = "encoding-type") <- "dict"
  hdf5r::h5attr(varp, which = "encoding-version") <- "0.1.0"

  # Close the .h5ad file and its groups
  adata_new$close_all()

  cli_alert_success("Exported anndata file")
}


#' Utility function used to merge and write the metadata.json file
#'
#' @param fs_map A \code{tbl_df} object with information about the component IDs
#' and what pixel file they belong to
#' @param pxl_folder A character string specifying the path to a temporary
#' directory where the merged metadata.json file will be created.
#'
#' @return Nothing
#'
#' @noRd
#'
.write_json_meta_data <- function (
  fs_map,
  pxl_folder
) {
  # If there is only one pixel file, we just copy the metadata.json file
  if (nrow(fs_map) == 1) {
    f <- fs_map$pxl_file
    unzip(f, exdir = pxl_folder, files = "metadata.json")
  } else {
    # Read json files from multiple pixel files and merge them into a single json file
    json_files <- sapply(seq_len(nrow(fs_map)), function(i) {
      f <- fs_map$pxl_file[i]
      unzip(f, exdir = pxl_folder, files = "metadata.json")
      sample_json <- file.path(pxl_folder, paste0("sample", i, ".json"))
      fs::file_move(file.path(pxl_folder, "metadata.json"), sample_json)
      return(sample_json)
    })
    merged_sample_json <- list(samples = lapply(json_files, function(f) {
      jsonlite::read_json(f)
    }) %>%
      set_names(nm = paste0("sample", seq_along(json_files))))
    jsonlite::write_json(merged_sample_json, file.path(pxl_folder, "metadata.json"), auto_unbox = TRUE)
    fs::file_delete(json_files)
  }
  cli_alert_success("Exported merged meta data")
}


#' Utility function used to write spatial metrics to a .parquet file
#'
#' @param object A \code{Seurat} object
#' @param pxl_folder A character string specifying the path to a temporary
#' directory where the polarization.parquet and colocalization.parquet files
#' will be created.
#'
#' @return Nothing
#'
#' @noRd
#'
.write_spatial_metrics_to_parquet <- function (
  object,
  pxl_folder
) {

  sb <- cli_status("{symbol$arrow_right} Writing spatial metrics to .parquet files...")
  cli_status_update(id = sb,
                    "{symbol$arrow_right} Writing polarization scores to polarization.parquet")

  # Fetch polarization scores from Seurat object and write to a .parquet file
  PolarizationScores(object) %>%
    arrow::write_dataset(hive_style = FALSE, format = "parquet",
                         compression = "zstd", path = pxl_folder,
                         max_rows_per_group = nrow(.))

  # The file is written as part-0.parquet, so we need to rename it to polarization.parquet
  fs::file_move(file.path(pxl_folder, "part-0.parquet"),
                file.path(pxl_folder, paste0("polarization.parquet")))


  cli_status_update(id = sb,
                    "{symbol$arrow_right} Writing colocalization scores to colocalization.parquet")

  # Fetch polarization scores from Seurat object and write to a .parquet file
  ColocalizationScores(object) %>%
    arrow::write_dataset(hive_style = FALSE, format = "parquet",
                         compression = "zstd", path = pxl_folder,
                         max_rows_per_group = nrow(.))

  # The file is written as part-0.parquet, so we need to rename it to polarization.parquet
  fs::file_move(file.path(pxl_folder, "part-0.parquet"),
                file.path(pxl_folder, paste0("colocalization.parquet")))

  cli_status_clear(id = sb)
  cli_alert_success("Exported spatial metrics")
}


#' Utility function used to merge and write the edgelists to a .parquet file
#'
#' @param fs_map A \code{tbl_df} object with information about the component IDs
#' and what pixel file they belong to
#' @param pxl_folder A character string specifying the path to a temporary
#' directory where the merged edgelist.parquet file will be created.
#'
#' @return Nothing
#'
#' @noRd
#'
.merge_and_write_edgelists_to_parquet <- function (
  fs_map,
  pxl_folder
) {

  # Unzip and collect edgelist files
  sb <- cli_status("{symbol$arrow_right} Collecting edge lists from {nrow(fs_map)} files")
  parquet_files <- sapply(seq_len(nrow(fs_map)), function(i) {
    f <- fs_map$pxl_file[i]
    unzip(f, exdir = pxl_folder, files = "edgelist.parquet")
    sample_parquet <- file.path(pxl_folder, paste0("sample", i, ".parquet"))
    fs::file_move(file.path(pxl_folder, "edgelist.parquet"), sample_parquet)
    cli_status_update(id = sb,
                      "{symbol$arrow_right} Extracting sample {i} edge list")
    return(sample_parquet)
  })

  # Unnset fs_map
  fs_map_unnest <- fs_map %>%
    tidyr::unnest(cols = "id_map")

  # Update component IDs in edgelist parquet files
  cli_status_update(id = sb,
                    "{symbol$arrow_right}  Merging edge lists")
  all_edgelists <- lapply(seq_along(parquet_files), function(i) {

    f <- parquet_files[i]

    # Read edgelist data from the current sampe
    data <- arrow::open_dataset(sources = f)
    data$metadata

    # Get type for component
    component_dtype <- arrow::schema(data)["component"]$component$type$name

    # Fetch the current_id and original_id to convert the component IDs
    # in the edgelist tables
    fs_map_unnest_filtered <- fs_map_unnest %>%
      filter(sample == i) %>%
      select(current_id, original_id) %>%
      rename(component = original_id) %>%
      {
        if (component_dtype == "large_utf8") {
          mutate(., across(current_id, factor))
        } else {
          mutate(., across(everything(), factor))
        }
      }

    # Convert fs_map_unnest_filtered to an arrow table
    a_table <- arrow::arrow_table(fs_map_unnest_filtered,
                                  schema = arrow::unify_schemas(
                                    arrow::schema(current_id = arrow::dictionary(arrow::int32(), arrow::utf8())),
                                    arrow::schema(data)["component"]
                                  ))

    # Update component IDs in the edgelist table using a join operation
    data <- data %>%
      left_join(a_table, by = "component") %>%
      select(-component) %>%
      rename(component = current_id) %>%
      # Filter out components that are not in the current_id
      filter(component %in% (fs_map_unnest_filtered$current_id %>% levels())) %>%
      # Force computation (this can be slow and requires the data to be loaded in memory)
      compute()

    cli_status_update(id = sb,
                      "{symbol$arrow_right} Loaded sample {i} edge list")

    return(data)
  })

  # Merge arrow tables
  cli_status_update(id = sb,
                    "{symbol$arrow_right} Merging edge lists")
  data_merged <- Reduce(arrow::concat_tables, all_edgelists)
  cli_status_update(id = sb,
                    "{symbol$arrow_right} Exporting merged edge list")

  # Write parquet file
  arrow::write_parquet(data_merged, file.path(pxl_folder, "edgelist.parquet"))
  cli_status_clear(id = sb)
  cli_alert_success("Exported merged edge list")

  # Remove sample parquet files
  fs::file_delete(parquet_files)
}


#' @param cg_list A list of \code{CellGraph} objects
#' @param pxl_folder A character string specifying the path to a temporary
#' directory where the merged layouts.parquet directory will be created.
#'
#' @return Nothing
#'
#' @noRd
#'
.merge_layout_with_counts_and_write_to_parquet <- function (
  cg_list,
  pxl_folder
) {

  sb <- cli_status("{symbol$arrow_right} Fetching layouts and marker count data for {length(cg_list)} components...")

  all_data <- lapply(names(cg_list), function(nm) {

    cg <- cg_list[[nm]]

    if (is.null(cg@layout)) return(invisible(NULL))
    if (is.null(cg@counts)) return(invisible(NULL))

    layouts <- cg@layout

    layout_types <- names(layouts)
    layout_types <- sapply(layout_types, function(layout_type) {
      ifelse((!str_detect(layout_type, "_3d$")) & (ncol(layouts[[layout_type]]) == 3),
             paste0(layout_type, "_3d"),
             layout_type)
    })
    layouts <- set_names(layouts, layout_types %>% unname())

    merged_counts <- lapply(layouts, function(layout) {

      # Combine counts with layout
      if (nrow(cg@counts) != nrow(layout))
        abort("Counts and layout must have the same number of rows")
      return(cg@counts %>% t())
    }) %>% do.call(cbind, .)
    #merged_counts <- SeuratObject::RowMergeSparseMatrices(all_counts[[1]], all_counts[-1])

    merged_layouts <- lapply(names(layouts), function(layout_type) {
      layout_tbl <- layouts[[layout_type]]
      if (ncol(layout_tbl) == 2 & !"z" %in% names(layout_tbl)) {
        layout_tbl <- layout_tbl %>% mutate(z = NA_real_)
      }
      layout_tbl <- layout_tbl[, c("x", "y", "z")]
      layout_tbl <- layout_tbl %>%
        mutate(layout = layout_type, component = nm, graph_projection = attr(cg@cellgraph, "type"))
      return(layout_tbl)
    }) %>% do.call(bind_rows, .)

    return(list(merged_counts = merged_counts, merged_layouts = merged_layouts))

  })

  cli_status_update(id = sb,
                    "{symbol$arrow_right} Merging marker count data")

  # merge all count data
  all_count_data_merged <- lapply(all_data, function(data) data$merged_counts)
  all_count_data_merged <- SeuratObject::RowMergeSparseMatrices(all_count_data_merged[[1]], all_count_data_merged[-1]) %>% t()

  cli_status_update(id = sb,
                    "{symbol$arrow_right} Merging layouts")

  # merge all layouts
  all_layout_merged <- lapply(all_data, function(data) data$merged_layouts) %>% do.call(bind_rows, .)

  cli_status_update(id = sb,
                    "{symbol$arrow_right} Merging count data with layouts")

  # Set schema for arrow table
  # Marker counts are "int64" and layout coordinates are "double"
  arr_table_schema <- arrow::schema(
    c(rep(list(arrow::int64()), ncol(all_count_data_merged)),
      c(rep(list(double()), ncol(all_layout_merged) - 3),
        arrow::string(),
        arrow::string(),
        arrow::string())) %>%
      set_names(nm = c(colnames(all_count_data_merged), colnames(all_layout_merged))))

  # Merge data
  all_layout_and_counts_merged <-
    cbind(all_count_data_merged %>% as.matrix(), all_layout_merged)

  # Convert data to arrow table
  layout_and_counts_arr <- arrow::arrow_table(all_layout_and_counts_merged,
                                             schema = arr_table_schema)

  # Group by graph_projection, layout and component
  layout_and_counts_arr <- layout_and_counts_arr %>%
    group_by(graph_projection, layout, component)

  # Write dataset
  arrow::write_dataset(layout_and_counts_arr,
                       path = file.path(pxl_folder, "layout.parquet"))

  cli_status_clear(id = sb)
  cli_alert_success("Exported layouts")
}
