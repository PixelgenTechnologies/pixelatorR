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


#' Fetches the config values to be used with DuckDB
#'
#' Configuration options can be found here:
#' https://duckdb.org/docs/stable/configuration/overview#global-configuration-options
#'
#' @return a named list with the config options to be used
#'
get_duckdb_config <- function() {
  config <- list()
  temp_dir <- Sys.getenv("PIXELATOR_DUCKDB_TEMP_DIR")
  if (temp_dir != "") {
    config$temp_directory <- temp_dir
  }
  config
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
              config = get_duckdb_config(),
              read_only = TRUE
            )
        },
        silent = TRUE
      )
      if (inherits(code, "try-error")) {
        cli::cli_abort(
          c(
            "i" = "The PXL file must contain PNA data.",
            "i" = "The PXL file might be blocked by another connection",
            "x" = "Failed to connect to the database."
          )
        )
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
              config = get_duckdb_config(),
              read_only = TRUE
            )
        },
        silent = TRUE
      )
      if (inherits(code, "try-error")) {
        cli::cli_abort(
          c(
            "i" = "The PXL file must contain PNA data.",
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
    #'   columns_filter = list("component" = c("0a45497c6bfbfb22", "2708240b908e2eba"))
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
        tibble::column_to_rownames("index") %>%
        as.matrix() %>%
        as("dgCMatrix") %>%
        Matrix::t()
    },
    #' @description
    #' Fetches the proximity scores from the database
    #'
    #' @param calc_log2_ratio A logical specifying whether to calculate and add
    #' a log2ratio column to the output table. Default is \code{TRUE}
    #' @param lazy A logical specifying whether to load the data lazily. If \code{TRUE},
    #' a \code{tbl_lazy} object is returned.
    #'
    #' @examples
    #' # Fetch the proximity scores
    #' prox <- db$proximity()
    #' prox %>% head()
    #'
    #' @return A \code{tbl_df} with the proximity scores
    #'
    proximity = function(calc_log2_ratio = TRUE, lazy = FALSE) {
      assert_single_value(calc_log2_ratio, "bool")
      assert_single_value(lazy, "bool")

      self$check_connection()
      if (lazy) {
        proximity_scores <- tbl(private$con, "proximity")
      } else {
        proximity_scores <- self$fetch_table("proximity")
      }
      proximity_scores <- proximity_scores %>%
        select(-matches("__index"))

      # Calculate log2ratio
      if (calc_log2_ratio) {
        if (!"join_count" %in% colnames(proximity_scores)) {
          cli::cli_abort(
            c("x" = "The Proximity score table does not contain the required {.str join_count} column")
          )
        }
        if (!"join_count_expected_mean" %in% colnames(proximity_scores)) {
          cli::cli_abort(
            c("x" = "The Proximity score table does not contain the required {.str join_count_expected_mean} column")
          )
        }
        proximity_scores <- proximity_scores %>%
          mutate(log2_ratio = log2(pmax(join_count, 1) / pmax(join_count_expected_mean, 1)))
      }

      proximity_scores <- proximity_scores %>% compute()

      if (!lazy) {
        proximity_scores <- proximity_scores %>% as_tibble()
      }

      return(proximity_scores)
    },
    #' @description
    #' Compute analytical proximity scores directly from the edgelist
    #'
    #' The calculation is performed entirely in DuckDB and does not require
    #' loading \code{CellGraph} objects into memory. Since PXL files are opened
    #' read-only, the result is materialized as a temporary table that exists
    #' for the lifetime of this connection. By default, the temporary table is
    #' named \code{proximity} and shadows a persisted table with the same name.
    #'
    #' By default, the method returns observed and expected join counts together
    #' with the \code{log2_ratio}. Z-scores and p-values are optional because they
    #' require DuckDB's stochastic extension.
    #'
    #' @param components A character vector of component names to include, or
    #' \code{NULL} to include all components.
    #' @param markers A character vector of marker names to include, or
    #' \code{NULL} to include all markers.
    #' @param calc_z_score Logical indicating whether to calculate z-scores and
    #' p-values. Default is \code{FALSE}.
    #' @param min_marker_count An integer specifying the minimum total marker count
    #'   (sum of UMI counts across umi1 and umi2) for a marker to be considered in proximity pair
    #'   expansion. Defaults to \code{10L}. Set to \code{NULL} or \code{0L} to include all markers.
    #'   Note that all edges are retained for computing expected frequencies and observed join counts
    #'   to prevent bias.
    #' @param batch_size An integer indicating the number of components to
    #' process in each batch, or \code{NULL} to process all requested components
    #' in a single query. Defaults to \code{100}.
    #' @param name Name of the temporary table in which to store the results.
    #'
    #' @examples
    #' db <- PixelDB$new(pxl_file)
    #' proximity <- db$compute_proximity_scores(
    #'   components = "0a45497c6bfbfb22",
    #'   markers = c("B2M", "HLA-ABC")
    #' )
    #' proximity %>% head()
    #'
    #' @return A \code{tbl_lazy} referring to the temporary proximity table
    #'
    compute_proximity_scores = function(
      components = NULL,
      markers = NULL,
      calc_z_score = FALSE,
      min_marker_count = 10L,
      batch_size = 100L,
      name = "proximity"
    ) {
      self$check_connection()
      assert_vector(components, "character", n = 1, allow_null = TRUE)
      assert_vector(markers, "character", n = 1, allow_null = TRUE)
      assert_single_value(calc_z_score, "bool")
      assert_single_value(min_marker_count, "int", allow_null = TRUE)
      if (!is.null(min_marker_count) && min_marker_count < 0L) {
        cli::cli_abort(c("x" = "{.arg min_marker_count} must be a non-negative integer or NULL."))
      }
      assert_single_value(batch_size, "int", allow_null = TRUE)
      if (!is.null(batch_size) && batch_size < 1L) {
        cli::cli_abort(c("x" = "{.arg batch_size} must be a positive integer."))
      }
      assert_single_value(name, "string")

      if (calc_z_score) {
        # The stochastic extension supplies dist_normal_cdf(), which is used for
        # the two-sided p-value. It may already be installed by DuckDB.
        tryCatch(
          DBI::dbExecute(private$con, "LOAD stochastic"),
          error = function(e) {
            DBI::dbExecute(private$con, "INSTALL stochastic FROM community")
            DBI::dbExecute(private$con, "LOAD stochastic")
          }
        )
      }

      edgelist_columns <- DBI::dbListFields(private$con, "edgelist")
      required_columns <- c("component", "umi1", "umi2", "marker_1", "marker_2")
      missing_columns <- setdiff(required_columns, edgelist_columns)
      if (length(missing_columns) > 0) {
        cli::cli_abort(
          c(
            "x" = "The edgelist is missing required column{?s}: {.field {missing_columns}}."
          )
        )
      }

      # If components is NULL and batching is enabled, resolve all component IDs
      all_components <- components
      if (!is.null(batch_size)) {
        if (is.null(all_components)) {
          all_tables <- DBI::dbListTables(private$con)
          if ("__adata__obs" %in% all_tables) {
            all_components <- DBI::dbGetQuery(
              private$con,
              "SELECT index AS component FROM __adata__obs"
            )[["component"]]
          } else {
            all_components <- DBI::dbGetQuery(
              private$con,
              "SELECT DISTINCT component FROM edgelist"
            )[["component"]]
          }
        }
      }

      # If batching is applicable (more components than batch_size)
      if (!is.null(batch_size) && length(all_components) > batch_size) {
        component_batches <- split(
          all_components,
          ceiling(seq_along(all_components) / batch_size)
        )
        n_batches <- length(component_batches)
        table_name <- as.character(DBI::dbQuoteIdentifier(private$con, name))

        pb <- cli_progress_bar(
          "Computing proximity scores",
          total = n_batches,
          .auto_close = FALSE
        )

        for (b in seq_along(component_batches)) {
          batch_comps <- component_batches[[b]]
          batch_query <- private$build_proximity_query(
            components = batch_comps,
            markers = markers,
            calc_z_score = calc_z_score,
            min_marker_count = min_marker_count,
            edgelist_columns = edgelist_columns
          )

          if (b == 1L) {
            DBI::dbExecute(
              private$con,
              glue::glue(
                "CREATE OR REPLACE TEMPORARY TABLE {table_name} AS ",
                "{batch_query}"
              )
            )
          } else {
            DBI::dbExecute(
              private$con,
              glue::glue(
                "INSERT INTO {table_name} ",
                "{batch_query}"
              )
            )
          }
          cli_progress_update(id = pb)
        }
        cli_progress_done(id = pb)

        return(tbl(private$con, name))
      }

      analysis_query <- private$build_proximity_query(
        components = components,
        markers = markers,
        calc_z_score = calc_z_score,
        min_marker_count = min_marker_count,
        edgelist_columns = edgelist_columns
      )

      table_name <- as.character(DBI::dbQuoteIdentifier(private$con, name))
      DBI::dbExecute(
        private$con,
        glue::glue(
          "CREATE OR REPLACE TEMPORARY TABLE {table_name} AS ",
          "{analysis_query}"
        )
      )

      return(tbl(private$con, name))
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
    #' @param lazy A logical specifying whether to load the data lazily. If \code{TRUE},
    #' a \code{tbl_lazy} object is returned.
    #' @param include_all_columns Logical specifying whether to include all columns in the output.
    #'
    #' @examples
    #' # Fetch edgelists
    #' db$components_edgelist("0a45497c6bfbfb22") %>% head()
    #'
    #' @return A \code{data.frame} with the component edgelist:
    #'  - umi1: A unique ID of the first RCA product
    #'  - umi2: A unique ID of the second RCA product
    #'  - marker_1: The first protein
    #'  - marker_2: The second protein
    #'  - component: The component name
    #'
    #' If \code{include_all_columns = TRUE}, the following columns are also returned:
    #'  - read_count: The number of reads supporting the edge
    #'  - uei_count: The number of unique event identifiers (UEIs) supporting the edge.
    #'    This column is optional and is only returned if it is present in the edgelist table.
    #'
    components_edgelist = function(components,
                                   umi_data_type = c("int64", "string", "suffixed_string"),
                                   lazy = FALSE,
                                   include_all_columns = FALSE) {
      self$check_connection()
      assert_vector(components, "character", n = 1, allow_null = TRUE)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)
      assert_single_value(lazy, type = "bool")
      assert_single_value(include_all_columns, type = "bool")

      umi_data_type <- match.arg(umi_data_type, choices = c("int64", "string", "suffixed_string"))

      el <- tbl(private$con, "edgelist")

      if (!is.null(components)) {
        el <- el %>%
          filter(component %in% components)
      }

      if (include_all_columns) {
        el <- el %>%
          mutate(across(any_of(c("read_count", "uei_count")), as.integer))
      } else {
        el <- el %>%
          select(-any_of(c("read_count", "uei_count")))
      }

      if (umi_data_type == "string") {
        el <- el %>%
          mutate(
            umi1 = as.character(umi1),
            umi2 = as.character(umi2)
          )
      }
      if (umi_data_type == "suffixed_string") {
        el <- el %>%
          mutate(
            umi1 = as.character(umi1) %>% stringr::str_c("-umi1"),
            umi2 = as.character(umi2) %>% stringr::str_c("-umi2")
          )
      }

      if (lazy) {
        # Register as a temporary VIEW instead of materializing a table
        sql_query <- dbplyr::sql_render(el)
        DBI::dbExecute(
          private$con,
          paste0("CREATE OR REPLACE TEMPORARY VIEW edgelist_modified AS ", sql_query)
        )

        el <- tbl(private$con, "edgelist_modified")
      } else {
        el <- el %>% collect()
      }
      return(el)
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
    #' db$components_layout("0a45497c6bfbfb22")[[1]] %>% head()
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
    #' db$components_marker_counts("0a45497c6bfbfb22")[[1]][1:3, 1:4]
    #'
    #' @return A list with \code{dgCMatrix} matrices or \code{tbl_df}'s with the marker counts
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
        data <- component_query(x) %>%
          DBI::dbGetQuery(private$con, .) %>%
          mutate(across(everything(), ~ factor(.x, levels = unique(.x))))
        m <- Matrix::sparseMatrix(
          i = as.integer(data$name),
          j = as.integer(data$marker),
          x = 1,
        )
        rownames(m) <- levels(data$name)
        colnames(m) <- levels(data$marker)
        return(m)
      }) %>%
        set_names(nm = components)

      if (!as_sparse) {
        component_list <- lapply(component_list, function(x) {
          x %>%
            as.matrix() %>%
            as_tibble(rownames = "name")
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
    export_parquet = function(parquet_file,
                              table_name = c("proximity", "edgelist", "layouts"),
                              compression = c("snappy", "zstd"),
                              compression_level = 1L) {
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
    con = NULL,
    build_proximity_query = function(components, markers, calc_z_score, min_marker_count, edgelist_columns) {
      group_columns <- c(
        if ("sample" %in% edgelist_columns) "sample",
        "component"
      )
      group_sql <- paste(group_columns, collapse = ", ")
      qualified_group_sql <- function(alias) {
        paste0(alias, ".", group_columns, collapse = ", ")
      }
      join_group_sql <- function(left, right) {
        paste0(left, ".", group_columns, " = ", right, ".", group_columns,
          collapse = " AND "
        )
      }
      am_group_sql <- qualified_group_sql("am")
      t1_group_sql <- qualified_group_sql("t1")
      am_s1_join_sql <- join_group_sql("am", "s1")
      am_s2_join_sql <- join_group_sql("am", "s2")
      t1_t2_join_sql <- join_group_sql("t1", "t2")
      t1_ge_join_sql <- join_group_sql("t1", "ge")

      sql_values <- function(values) {
        paste(DBI::dbQuoteString(private$con, unique(values)), collapse = ", ")
      }
      component_filter <- if (is.null(components)) {
        "TRUE"
      } else {
        glue::glue("component IN ({sql_values(components)})")
      }
      marker_filter_expected <- if (is.null(markers)) {
        "TRUE"
      } else {
        glue::glue("marker IN ({sql_values(markers)})")
      }
      marker_filter_observed <- if (is.null(markers)) {
        "TRUE"
      } else {
        marker_values <- sql_values(markers)
        glue::glue(
          "marker_1 IN ({marker_values}) ",
          "AND marker_2 IN ({marker_values})"
        )
      }

      min_marker_filter_sql <- if (!is.null(min_marker_count) && min_marker_count > 0L) {
        glue::glue(
          "WHERE (t1.marker_1_count + t1.marker_2_count) >= {min_marker_count} ",
          "AND (t2.marker_1_count + t2.marker_2_count) >= {min_marker_count}"
        )
      } else {
        ""
      }

      coalesced_group_sql <- paste0(
        "exp.", group_columns, " AS ",
        group_columns,
        collapse = ", "
      )
      final_join_sql <- join_group_sql("exp", "obs")
      if (calc_z_score) {
        expected_sd_sql <- paste(
          ",",
          "CASE",
          "  WHEN t1.marker = t2.marker THEN",
          "    SQRT(t1.f_umi1 * t2.f_umi2 * (1.0 - (t1.f_umi1 * t2.f_umi2)) * ge.n_edges)",
          "  ELSE",
          "    SQRT(",
          "      (t1.f_umi1 * t2.f_umi2 * (1.0 - (t1.f_umi1 * t2.f_umi2)) +",
          "       t2.f_umi1 * t1.f_umi2 * (1.0 - (t2.f_umi1 * t1.f_umi2))) * ge.n_edges",
          "    )",
          "END AS join_count_expected_sd"
        )
        result_sd_sql <- paste(
          ",",
          "GREATEST(COALESCE(exp.join_count_expected_sd, 0), 1e-6)",
          "AS join_count_expected_sd"
        )
        z_score_sql <- paste(
          "(join_count - join_count_expected_mean) / join_count_expected_sd",
          "AS join_count_z,",
          "2 * (",
          "  1 - dist_normal_cdf(",
          "    0.0,",
          "    join_count_expected_sd,",
          "    ABS(join_count - join_count_expected_mean)",
          "  )",
          ") AS join_count_p,"
        )
      } else {
        expected_sd_sql <- ""
        result_sd_sql <- ""
        z_score_sql <- ""
      }

      glue::glue(
        "
        WITH
        current_edgelist AS (
          SELECT {group_sql}, umi1, umi2, marker_1, marker_2
          FROM edgelist
          WHERE {component_filter}
        ),
        group_edges AS (
          SELECT {group_sql}, COUNT(*) AS n_edges
          FROM current_edgelist
          GROUP BY {group_sql}
        ),
        unique_m1 AS (
          SELECT DISTINCT {group_sql}, umi1, marker_1 FROM current_edgelist
        ),
        stats_m1 AS (
          SELECT
            {group_sql},
            marker_1 AS marker,
            COUNT(*) AS marker_1_count,
            COUNT(*) * 1.0 /
              SUM(COUNT(*)) OVER (PARTITION BY {group_sql}) AS f_umi1
          FROM unique_m1
          GROUP BY {group_sql}, marker_1
        ),
        unique_m2 AS (
          SELECT DISTINCT {group_sql}, umi2, marker_2 FROM current_edgelist
        ),
        stats_m2 AS (
          SELECT
            {group_sql},
            marker_2 AS marker,
            COUNT(*) AS marker_2_count,
            COUNT(*) * 1.0 /
              SUM(COUNT(*)) OVER (PARTITION BY {group_sql}) AS f_umi2
          FROM unique_m2
          GROUP BY {group_sql}, marker_2
        ),
        all_markers AS (
          SELECT {group_sql}, marker FROM stats_m1 WHERE {marker_filter_expected}
          UNION
          SELECT {group_sql}, marker FROM stats_m2 WHERE {marker_filter_expected}
        ),
        combined_stats AS (
          SELECT
            {am_group_sql},
            am.marker,
            COALESCE(s1.f_umi1, 0.0) AS f_umi1,
            COALESCE(s2.f_umi2, 0.0) AS f_umi2,
            COALESCE(s1.marker_1_count, 0) AS marker_1_count,
            COALESCE(s2.marker_2_count, 0) AS marker_2_count
          FROM all_markers am
          LEFT JOIN stats_m1 s1
            ON {am_s1_join_sql}
            AND am.marker = s1.marker
          LEFT JOIN stats_m2 s2
            ON {am_s2_join_sql}
            AND am.marker = s2.marker
        ),
        observed_agg AS (
          SELECT
            {group_sql},
            LEAST(marker_1, marker_2) AS marker_A,
            GREATEST(marker_1, marker_2) AS marker_B,
            COUNT(*) AS join_count
          FROM current_edgelist
          WHERE {marker_filter_observed}
          GROUP BY {group_sql}, marker_A, marker_B
        ),
        expected_agg AS (
          SELECT
            {t1_group_sql},
            t1.marker AS marker_A,
            t2.marker AS marker_B,
            CASE
              WHEN t1.marker = t2.marker THEN t1.f_umi1 * t2.f_umi2 * ge.n_edges
              ELSE (t1.f_umi1 * t2.f_umi2 + t2.f_umi1 * t1.f_umi2) * ge.n_edges
            END AS join_count_expected_mean
            {expected_sd_sql}
          FROM combined_stats t1
          JOIN combined_stats t2
            ON {t1_t2_join_sql}
            AND t1.marker <= t2.marker
          JOIN group_edges ge
            ON {t1_ge_join_sql}
          {min_marker_filter_sql}
        ),
        results AS (
          SELECT
            {coalesced_group_sql},
            exp.marker_A AS marker_1,
            exp.marker_B AS marker_2,
            COALESCE(obs.join_count, 0) AS join_count,
            exp.join_count_expected_mean
            {result_sd_sql}
          FROM expected_agg exp
          LEFT JOIN observed_agg obs
            ON {final_join_sql}
            AND obs.marker_A = exp.marker_A
            AND obs.marker_B = exp.marker_B
        )
        SELECT
          *,
          {z_score_sql}
          LOG2(
            GREATEST(join_count, 1) /
            GREATEST(join_count_expected_mean, 1)
          ) AS log2_ratio
        FROM results
        "
      )
    },
    finalize = function() {
      if (DBI::dbIsValid(private$con)) {
        DBI::dbDisconnect(private$con)
      }
    }
  ),
  cloneable = FALSE
)
