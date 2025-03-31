#' Lazy load merged tables from multiple databases (PXL files)
#'
#' @param fs_map A \code{tbl_df} with information about the PXL file paths,
#' and component IDs.
#' @param table_name The name of the table to load, e.g. "proximity".
#' @param calc_log2ratio A logical indicating whether to calculate the
#' log2 ratio proximity score. This is only used for the "proximity" table.
#' @param call The calling environment.
#'
#' @return A \code{tbl_lazy} table.
#'
#' @examples
#' # SETUP
#' tmp_dir <- file.path(fs::path_temp(), "test")
#' fs::dir_create(tmp_dir)
#' for (i in 1:3) {
#'   fs::file_copy(
#'     minimal_pna_pxl_file(),
#'     file.path(tmp_dir, glue::glue("S{i}.pxl")),
#'     overwrite = TRUE
#'   )
#' }
#' pxl_files <- list.files(tmp_dir, full.names = TRUE)
#'
#' # Create a merged Seurat object
#' se_list <- lapply(seq_along(pxl_files), function(i) {
#'   se <- ReadPNA_Seurat(pxl_files[i])
#'   se$sample_id <- glue::glue("Sample{i}")
#'   return(se)
#' })
#' se <- merge(se_list[[1]], se_list[-1], add.cell.ids = LETTERS[1:3])
#' # SETUP END
#'
#' # Lazy load merged proximity tables
#' proximity_lazy <- pixelatorR:::.lazy_load_table(FSMap(se), "proximity")
#' proximity_lazy
#'
#' # Lazy load merged edgelists
#' edgelist_lazy <- pixelatorR:::.lazy_load_table(FSMap(se), "edgelist")
#' edgelist_lazy
#'
#' \dontrun{
#' # Lazy load merged layotus
#' edgelist_lazy <- pixelatorR:::.lazy_load_table(FSMap(se), "layouts")
#' edgelist_lazy
#' }
#'
#' @noRd
#'
.lazy_load_table <- function(
  fs_map,
  table_name = "proximity",
  calc_log2ratio = TRUE,
  call = caller_env()
) {
  # Check that PXL files are unique
  dup_pxl_files <- fs_map$pxl_file[duplicated(fs_map$pxl_file)]
  if (length(dup_pxl_files) > 0) {
    cli::cli_abort(
      c(
        "i" = "All PXL files must be unique.",
        "x" = "The following PXL file(s) are duplicated: ",
        "x" = "{.file {dup_pxl_files}}"
      )
    )
  }

  # Create an in-memory federated database to create a view of the combined tables
  con_federated <- DBI::dbConnect(duckdb::duckdb(bigint = "integer64"), ":memory:")

  # Attach the source databases
  for (i in seq_len(nrow(fs_map))) {
    # Attach each PXL file database
    DBI::dbExecute(con_federated, glue::glue("ATTACH '{fs_map$pxl_file[i]}' AS db{i} (TYPE duckdb, READ_ONLY)"))
    # Copy the ID map to the federated database
    copy_to(con_federated, fs_map$id_map[[i]] %>% rename(component = original_id), glue::glue("id_map{i}"))
  }

  # Select column names to keep
  column_names <- DBI::dbGetQuery(con_federated, glue::glue("SELECT * FROM db1.{table_name} LIMIT 1")) %>% names()
  column_names <- column_names[!stringr::str_detect(column_names, "__index_level_0__")]
  column_names <- paste(setdiff(column_names, "component"), collapse = ", ")

  # Build SQL query to create a view of the combined tables
  dbs <- seq_len(nrow(fs_map))
  if (calc_log2ratio && table_name == "proximity") {
    select_sql <- glue::glue(
      "SELECT {column_names}, ",
      "log2(GREATEST(join_count, 1.0)/GREATEST(join_count_expected_mean, 1.0)) ",
      "AS log2_ratio, "
    )
  } else {
    select_sql <- "SELECT {column_names}, "
  }
  union_sql <- paste(
    glue::glue(
      select_sql,
      "current_id AS component\n",
      "FROM db{dbs}.{table_name}\n",
      # Update component IDs using the ID map
      "JOIN id_map{dbs}\n",
      "ON db{dbs}.{table_name}.component = id_map{dbs}.component\n"
    ),
    # Combine all tables with UNION ALL
    collapse = "\nUNION ALL\n"
  )
  create_view_sql <- glue::glue("CREATE OR REPLACE VIEW combined_{table_name} AS {union_sql[1]}")

  # Execute the SQL query to create the view
  DBI::dbExecute(con_federated, create_view_sql)

  # Create a lazy table with dbplyr
  lazy_table <- tbl(con_federated, glue::glue("combined_{table_name}"))

  return(lazy_table)
}
