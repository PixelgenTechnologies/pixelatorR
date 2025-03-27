#' Assert that a file is a PXL file
#'
#' @param pxl_file A path to a PXL file
#'
#' @return Nothing
#'
#' @noRd
#'
assert_pxl_file <- function(pxl_file) {
  if (!fs::file_exists(pxl_file)) {
    cli::cli_abort(
      c(
        "x" = "File {.file {pxl_file}} does not exist"
      )
    )
  }
  if (!fs::path_ext(pxl_file) == "pxl") {
    cli::cli_abort(
      c(
        "i" = "{.var pxl_file} must be a path to a PXL file with extenstion {.val pxl}",
        "x" = "File {.file {pxl_file}} has extension {.val {fs::path_ext(pxl_file)}}"
      )
    )
  }
}


#' PXL database class
#'
#' This class provides an interface for working with a PXL file.
#' A PXL file is a duckdb database that contains various tables created
#' by the Pixelator data processing pipeline.
#'
#' The class provides methods to query the database and extract information
#' from it. The class is intended to be used internally by the package but
#' is exposed to the user for advanced use cases.
#'
#' Example usage:
#' ```r
#' # Create a PixelDB object with a connection to a PXL file
#' pxl_db <- PixelDB$new("path/to/pxl/file")
#' pxl_db$info()
#'
#' # Fetch count matrix
#' counts <- pxl_db$counts()
#'
#' # Fetch proximity scores
#' proximity <- pxl_db$proximity()
#' ```
#'
#' @export
#'
PixelDB <- R6Class(
  "PixelDB",
  public = list(
    #' @description
    #' Set up a connection to a PXL file
    #'
    #' @param file A path to a PXL file
    #'
    #' @examples
    #' library(dplyr)
    #'
    #' pxl_file <- minimal_pna_pxl_file()
    #' db <- PixelDB$new(pxl_file)
    #'
    #' @return An \code{R6} object representing the PXL database
    #'
    initialize = function(file) {
      private$file <- normalizePath(file)
      assert_pxl_file(file)
      code <- try(
        {
          private$con <-
            DBI::dbConnect(
              duckdb::duckdb(),
              bigint = "integer64",
              dbdir = private$file,
              read_only = TRUE
            )
        },
        silent = TRUE
      )
      if (inherits(code, "try-error")) {
        cli::cli_abort(
          c(
            "i" = "The PXL file must have been created with Pixelator version x.x.x or later.",
            "i" = "The PXL file might be blocked by another connection",
            "x" = "Failed to connect to the database."
          )
        )
      }
    },
    #' @description
    #' Cleanup method that is called when the \code{PixelDB} object is garbage collected
    #'
    #' @examples
    #' con <- db$.__enclos_env__$private$con
    #' rm(db)
    #' gc(full = FALSE)
    #'
    #' # Connection should now be closed
    #' DBI::dbIsValid(con) # FALSE
    #'
    #' @return Nothing
    #'
    finalize = function() {
      if (DBI::dbIsValid(private$con)) {
        DBI::dbDisconnect(private$con)
      }
    },
    #' @description
    #' Show information about tables in the PXL file
    #'
    #' @examples
    #' db <- PixelDB$new(pxl_file)
    #' db$info()
    #'
    #' @return A \code{tbl_df} with information about tables in the PXL file
    #'
    info = function() {
      DBI::dbGetQuery(private$con, "SHOW ALL TABLES") %>%
        as_tibble()
    },
    #' @description
    #' Send a database query
    #'
    #' @param sql An SQL query
    #'
    #' @examples
    #' # Select the proximity score table
    #' db$query("SELECT * FROM proximity") %>% head()
    #'
    #' @return A \code{data.frame} with the results of the query
    #'
    query = function(sql) {
      result <- try(
        {
          DBI::dbGetQuery(private$con, sql)
        },
        silent = TRUE
      )
      if (inherits(result, "try-error")) {
        cli::cli_abort(
          c(
            "x" = "Invalid query"
          )
        )
      }
      return(result)
    },
    #' @description
    #' Method to reconnect to the PXL file. This is useful if the connection is lost,
    #' for instance if the R session restarted.
    #'
    #' If the method fails, an error is thrown.
    #'
    #' @examples
    #' db$close()
    #' db$reconnect()
    #'
    #' @return Nothing
    #'
    reconnect = function() {
      if (DBI::dbIsValid(private$con)) {
        DBI::dbDisconnect(private$con) # Close old connection if needed
      }
      assert_pxl_file(private$file)
      code <- try(
        {
          private$con <-
            DBI::dbConnect(
              duckdb::duckdb(),
              bigint = "integer64",
              dbdir = private$file,
              read_only = TRUE
            )
        },
        silent = TRUE
      )
      if (inherits(code, "try-error")) {
        cli::cli_abort(
          c(
            "i" = "The PXL file must have been created with Pixelator version x.x.x or later.",
            "i" = "The PXL file might be blocked by another connection",
            "x" = "Failed to connect to the database."
          )
        )
      }
      cli::cli_inform(
        c("v" = "Connected")
      )
    },
    #' @description
    #' Check the connection to the PXL database
    #'
    #' If the connection is invalid, the \code{$reconnect()} method is called.
    #'
    #' @examples
    #' # If the connection is closed, the method will attempt to reconnect
    #' db$check_connection()
    #'
    #' @return Nothing
    #'
    check_connection = function() {
      if (is.null(private$con) || !DBI::dbIsValid(private$con)) {
        cli::cli_inform(
          c("i" = "Reconnecting")
        )
        self$reconnect()
      }
    },
    #' @description
    #' Get the names of the tables in the database
    #'
    #' @examples
    #' # Get the table names
    #' db$names()
    #'
    #' @return A character vector with table names that can be used
    #' in \code{$fetch_table()}
    #'
    names = function() {
      self$check_connection()
      self$info() %>%
        pull(name) %>%
        as.character()
    },
    #' @description
    #' Fetches an entire table from the database
    #'
    #' This general method can be used to fetch any table from the database.
    #' Large tables such as the edgelist are not recommended to be fetched in this way
    #' as you will get the entire table in memory.
    #'
    #' @param name The name of the table
    #'
    #' @examples
    #' # Fetch any table from the database
    #' db$fetch_table("proximity") %>% head()
    #'
    #' @return A \code{data.frame} with the table contents
    #'
    fetch_table = function(name) {
      self$check_connection()
      assert_single_value(name, "string")
      if (!name %in% self$names()) {
        cli::cli_abort(
          c(
            "i" = "{.arg name} should be one of {.val {self$names()}}",
            "x" = "You have attempted to fetch {.val {name}} which is not available"
          )
        )
      }
      DBI::dbGetQuery(
        private$con, glue("SELECT * FROM {name}")
      )
    },
    #' @description
    #' Fetches a subset of columns from a table in the database
    #'
    #' This method can be used to fetch a subset of columns from a table in the database.
    #'
    #' @param name The name of the table
    #' @param columns_filter A named list
    #'
    #' @examples
    #' # Fetch any table from the database and filter on the fly
    #' db$fetch_table_subset(
    #'   "proximity",
    #'   columns_filter = list("component" = c("3898b03349c6e28d", "406b9e5d80941ca0"))
    #' ) %>% head()
    #'
    #' @return A \code{data.frame} with the table contents
    #'
    fetch_table_subset = function(name, columns_filter) {
      self$check_connection()
      assert_single_value(name, "string")
      assert_class(columns_filter, "list")
      if (is.null(names(columns_filter))) {
        cli::cli_abort(
          c("x" = "{.arg columns_filter} must be a named {.cls list}")
        )
      }
      for (nm in names(columns_filter)) {
        if (!is.character(columns_filter[[nm]])) {
          cli::cli_abort(
            c(
              "i" = "Elements in {.arg columns_filter} must be {.cls character}",
              "x" = "Element {.val {nm}} in {.arg columns_filter} is {.cls {class(columns_filter[[nm]])}}"
            )
          )
        }
      }

      column_filter_collapsed <- sapply(names(columns_filter), function(nm) {
        x <- columns_filter[[nm]]
        glue("{nm} IN ({glue_sql('{x*}', .con = private$con)})")
      }) %>% paste(collapse = " AND ")
      sql_query <- glue("SELECT * FROM {name} WHERE {column_filter_collapsed}")
      DBI::dbGetQuery(private$con, sql_query)
    },
    #' @description
    #' Fetches the __adata__X count table from the database and converts it to a sparse \code{dgCMatrix}
    #'
    #' @examples
    #' # Fetch the antibody counts
    #' X <- db$counts()
    #' X[1:4, 1:4]
    #'
    #' @return A sparse matrix of class \code{dgCMatrix} with antibody counts
    #'
    counts = function() {
      self$check_connection()
      self$fetch_table("__adata__X") %>%
        mutate(across(where(bit64::is.integer64), as.integer)) %>%
        data.frame(row.names = 1, check.names = FALSE) %>%
        as.matrix() %>%
        as("dgCMatrix") %>%
        Matrix::t()
    },
    #' @description
    #' Fetches the proximity scores from the database
    #'
    #' @param calc_log2_ratio A logical specifying whether to calculate and add
    #' a log2ratio column to the output table. Default is \code{TRUE}
    #'
    #' @examples
    #' # Fetch the proximity scores
    #' prox <- db$proximity()
    #' prox %>% head()
    #'
    #' @return A \code{tbl_df} with the proximity scores
    #'
    proximity = function(calc_log2_ratio = TRUE) {
      self$check_connection()
      proximity_scores <- self$fetch_table("proximity") %>%
        as_tibble() %>%
        select(-matches("__index"))

      # Calculate log2ratio
      if (calc_log2_ratio) {
        if (!"join_count" %in% names(proximity_scores)) {
          cli::cli_abort(
            c("x" = "The Proximity score table does not contain the required {.str join_count} column")
          )
        }
        if (!"join_count_expected_mean" %in% names(proximity_scores)) {
          cli::cli_abort(
            c("x" = "The Proximity score table does not contain the required {.str join_count_expected_mean} column")
          )
        }
        proximity_scores <- proximity_scores %>%
          mutate(log2_ratio = log2(pmax(join_count, 1) / pmax(join_count_expected_mean, 1)))
      }

      return(proximity_scores)
    },
    #' @description
    #' Fetches the __adata__obs meta data
    #'
    #' @examples
    #' # Fetch the cell meta data
    #' db$cell_meta() %>% head()
    #'
    #' @return A \code{data.frame} with the cell meta data
    #'
    cell_meta = function() {
      self$check_connection()
      self$fetch_table("__adata__obs") %>%
        data.frame(row.names = 1, check.names = FALSE) %>%
        mutate(across(where(bit64::is.integer64), ~ as.integer(.x)))
    },
    #' @description
    #' Fetches the __adata__var meta data
    #'
    #' @examples
    #' # Fetch the protein meta data
    #' db$protein_meta() %>% head()
    #'
    #' @return A \code{data.frame} with the protein meta data
    #'
    protein_meta = function() {
      self$check_connection()
      self$fetch_table("__adata__var") %>%
        data.frame(row.names = 1, check.names = FALSE) %>%
        mutate(across(where(bit64::is.integer64), ~ as.integer(.x)))
    },
    #' @description
    #' Fetch the run meta data
    #'
    #' @examples
    #' # Fetch the run meta data
    #' db$run_meta()
    #'
    #' @return A \code{tbl_df} with the run meta data
    #'
    run_meta = function() {
      self$check_connection()
      string <- self$fetch_table("metadata")$value
      jsonlite::fromJSON(string) %>%
        as_tibble()
    },
    #' @description
    #' Fetches a component edgelist or the entire edgelist from the database
    #'
    #' The UMIs are encoded as int64 but since R doesn't support int64, the UMIs can be converted to
    #' character vectors be specifying the \code{umi_data_type}.
    #'
    #' @param components A character vector with component names or NULL to get all components
    #' @param umi_data_type One of "int64", "string" or "suffixed_string". Default is "int64".
    #'  - "int64": The UMIs are encoded as int64
    #'  - "string": The UMIs are encoded as character
    #'  - "suffixed_string": The UMIs are encoded as character with a suffix '-umi1' or '-umi2' added
    #' @param include_all_columns Logical specifying whether to include all columns in the output.
    #'
    #' @examples
    #' # Fetch edgelists
    #' db$components_edgelist("3898b03349c6e28d") %>% head()
    #'
    #' @return A \code{data.frame} with the component edgelist:
    #'  - umi1: A unique ID of the first RCA product
    #'  - umi2: A unique ID of the second RCA product
    #'  - marker_1: The first protein
    #'  - marker_2: The second protein
    #'  - component: The component name
    #'  - read_count: The number of reads supporting the edge
    #'  - uei_count: The number of unique event identifiers (UEIs) supporting the edge
    #'
    components_edgelist = function(
      components,
      umi_data_type = c("int64", "string", "suffixed_string"),
      include_all_columns = FALSE
    ) {
      self$check_connection()
      assert_vector(components, "character", n = 1, allow_null = TRUE)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)
      assert_single_value(include_all_columns, type = "bool")

      umi_data_type <- match.arg(umi_data_type, choices = c("int64", "string", "suffixed_string"))

      columns <- "marker_1, marker_2, component"
      if (include_all_columns) {
        columns <- glue::glue(
          "{columns}, ",
          "CAST(read_count AS SMALLINT) AS read_count, ",
          "CAST(uei_count AS SMALLINT) AS uei_count "
        )
      } else {
        columns <- glue::glue("{columns} ")
      }

      if (!is.null(components)) {
        components_sql <- glue_sql("{components*}", .con = private$con)
        components_filter_sql <- glue(" WHERE component IN ({components_sql})")
      } else {
        components_filter_sql <- ""
      }

      sql_query <- glue::glue(
        "SELECT ",
        switch(umi_data_type,
          "int64" = "umi1, ",
          "string" = "CAST(umi1 AS STRING) AS umi1, ",
          "suffixed_string" = "CAST(umi1 AS STRING) || '-umi1' AS umi1, "
        ),
        switch(umi_data_type,
          "int64" = "umi2, ",
          "string" = "CAST(umi2 AS STRING) AS umi2, ",
          "suffixed_string" = "CAST(umi2 AS STRING) || '-umi2' AS umi2, "
        ),
        columns,
        "FROM edgelist",
        components_filter_sql
      )
      DBI::dbGetQuery(private$con, sql_query)
    },
    #' @description
    #' Fetches layout for selected components
    #'
    #' This method fetches the x, y, z layout coordinates for selected PNA \code{components}.
    #'
    #' @param components A vector of component names or NULL for all components
    #' @param add_marker_counts Add marker counts in wide format to the layout tables
    #' @param verbose Print messages
    #'
    #' @examples
    #' # Fetch layouts (NOTE: This will only work if the layouts exist in the database)
    #' \dontrun{
    #' db$components_layout("3898b03349c6e28d")[[1]] %>% head()
    #' }
    #'
    #' @return A list with \code{tbl_df}'s with the layout coordinates and optionally marker counts
    #'
    components_layout = function(components, add_marker_counts = FALSE, verbose = TRUE) {
      self$check_connection()
      assert_vector(components, "character", allow_null = TRUE, n = 1)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)
      assert_single_value(add_marker_counts, "bool")

      # Construct the SQL query
      # The pixel_type column is used to determine the UMI suffix
      sql_query <- glue::glue(
        "SELECT CAST(index AS STRING) || ",
        "CASE WHEN pixel_type IN ('A', 'umi1') ",
        "THEN '-umi1' ELSE '-umi2' END AS ",
        "name, x, y, z FROM layouts"
      )

      components <- components %||% (self$counts() %>% colnames())

      if (verbose) cli::cli_alert_info("Fetching {.val {length(components)}} component layouts...")

      lapply_func <- ifelse(verbose, pbapply::pblapply, lapply)

      layout_list <- lapply_func(components, function(component) {
        sql_query <- glue::glue(
          sql_query,
          glue::glue(" WHERE component IN ('{component}')")
        )
        DBI::dbGetQuery(private$con, sql_query) %>%
          as_tibble()
      }) %>%
        set_names(nm = components)

      if (add_marker_counts) {
        if (verbose) cli::cli_alert_info("Fetching marker counts...")
        marker_counts <- self$components_marker_counts(components)
        if (verbose) cli::cli_alert_info("Adding marker counts to layout tables...")
        layout_list <- lapply_func(names(layout_list), function(nm) {
          layout_list[[nm]] %>%
            left_join(marker_counts[[nm]], by = c("name" = "name"))
        }) %>%
          set_names(nm = names(layout_list))
      }

      return(layout_list)
    },
    #' @description
    #' Fetches marker counts for components
    #'
    #' This method fetches the node marker counts for selected \code{components}. The node IDs
    #' are stored in the \code{name} column. If \code{components} is \code{NULL}, all
    #' component marker counts are returned with an additional \code{components} column in
    #' the resulting table.
    #'
    #' @param components A vector of component IDs or NULL for all components
    #' @param as_sparse Return the marker counts as a sparse matrix
    #' @param verbose Print messages
    #'
    #' @examples
    #' # Fetch marker counts
    #' db$components_marker_counts("3898b03349c6e28d")[[1]][1:3, 1:4]
    #'
    #' @return A list with \code{tbl_df}'s with the marker counts
    #'
    components_marker_counts = function(components, as_sparse = FALSE, verbose = FALSE) {
      self$check_connection()
      assert_vector(components, "character", allow_null = TRUE, n = 1)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)

      # Build the query for a single component
      component_query <- function(component) {
        component_filter_sql <- glue::glue("WHERE component IN ('{component}') ")
        sql_query <- glue::glue(
          "SELECT umi1 || '-umi1' as name, marker_1 as marker ", # Add -umi1 suffix
          "FROM edgelist ",
          component_filter_sql,
          "UNION ",
          "SELECT umi2 || '-umi2' as name, marker_2 as marker ", # Add -umi2 suffix
          "FROM edgelist ",
          component_filter_sql
        )
        return(sql_query)
      }

      components <- components %||% (self$counts() %>% colnames())

      lapply_func <- ifelse(verbose, pbapply::pblapply, lapply)

      component_list <- lapply_func(components, function(x) {
        component_query(x) %>%
          DBI::dbGetQuery(private$con, .) %>%
          as_tibble() %>%
          mutate(across(where(bit64::is.integer64), ~ as.integer(.x))) %>%
          mutate(n = 1) %>%
          pivot_wider(names_from = marker, values_from = n, values_fill = 0)
      }) %>%
        set_names(nm = components)

      if (as_sparse) {
        component_list <- lapply(component_list, function(x) {
          x %>%
            data.frame(row.names = 1, check.names = FALSE) %>%
            as.matrix() %>%
            as("dgCMatrix")
        })
      }

      return(component_list)
    },

    #' @description
    #' Export a table to a parquet file
    #'
    #' @param parquet_file Path to the parquet file
    #' @param table_name The name of the table to export
    #' @param compression The compression algorithm to use. Options are 'snappy' and 'zstd'.
    #' @param compression_level The compression level to use. Default is 1. Only used when
    #'  compression is 'zstd'.
    #'
    #' @examples
    #' # Export a table to a parquet file
    #' tmp_parquet_file <- fs::file_temp(ext = "parquet")
    #' db$export_parquet(tmp_parquet_file, "proximity")
    #' fs::file_exists(tmp_parquet_file)
    #'
    #' @return Nothing
    #'
    export_parquet = function(
      parquet_file,
      table_name = c("proximity", "edgelist", "layouts"),
      compression = c("snappy", "zstd"),
      compression_level = 1L
    ) {
      self$check_connection()
      assert_single_value(parquet_file, "string")
      assert_file_ext(parquet_file, "parquet")
      table_name <- match.arg(table_name, choices = c("proximity", "edgelist", "layouts"))
      compression <- match.arg(compression, c("snappy", "zstd"))
      assert_single_value(compression_level, "integer")
      assert_x_in_y(table_name, self$names())

      parquet_file <- normalizePath(parquet_file, mustWork = FALSE)

      if (compression != "snappy") {
        compression_sql <-
          glue::glue(
            "(FORMAT parquet, ",
            "COMPRESSION {compression}, ",
            "COMPRESSION_LEVEL {compression_level})"
          )
      } else {
        compression_sql <- "(FORMAT parquet)"
      }

      DBI::dbExecute(
        private$con,
        glue::glue(
          "COPY {table_name} ",
          "TO '{parquet_file}' ",
          compression_sql
        )
      )

      return(invisible(NULL))
    },
    #'
    #' @description
    #' Close connection
    #'
    #' @examples
    #' # Close the connection when finished
    #' db$close()
    #'
    #' @return Nothing
    #'
    close = function() {
      if (DBI::dbIsValid(private$con)) {
        DBI::dbDisconnect(private$con)
      }
    }
  ),
  private = list(
    file = NULL,
    con = NULL
  ),
  cloneable = FALSE
)
