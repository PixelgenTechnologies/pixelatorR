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
      private$component_index_cache <- list()
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

      # Check if the PXL file is marked as empty/null
      tables <- DBI::dbListTables(private$con)
      if ("metadata" %in% tables) {
        meta_query <- try(DBI::dbGetQuery(private$con, "SELECT value FROM metadata"), silent = TRUE)
        if (!inherits(meta_query, "try-error") && nrow(meta_query) > 0 && nzchar(meta_query$value[1])) {
          meta_parsed <- try(jsonlite::fromJSON(meta_query$value[1]), silent = TRUE)
          if (!inherits(meta_parsed, "try-error") && isTRUE(meta_parsed[["null"]])) {
            reason <- meta_parsed[["null_reason"]]
            self$close()
            # Leave `{reason}` for cli to substitute. Glue-formatting it first
            # would make cli parse braces inside the reason and throw a
            # different error than pixeldb_null_error.
            msg <- c("x" = "The PXL file is empty (marked as null).")
            if (is.character(reason) && length(reason) == 1L && !is.na(reason) && nzchar(reason)) {
              msg <- c(msg, "i" = "Reason: {reason}")
            }
            cli::cli_abort(msg, class = "pixeldb_null_error")
          }
        }
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
        private$con,
        glue("SELECT * FROM {name}")
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
      }) %>%
        paste(collapse = " AND ")
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
      if (!is.null(private$counts_cache)) {
        return(private$counts_cache)
      }
      df <- DBI::dbGetQuery(private$con, "SELECT * FROM __adata__X")
      rn <- df[["index"]]
      df[["index"]] <- NULL
      mat <- matrix(
        vapply(df, as.integer, integer(nrow(df))),
        nrow = nrow(df),
        ncol = length(df),
        dimnames = list(rn, names(df))
      )
      private$counts_cache <- Matrix::t(methods::as(mat, "dgCMatrix"))
      private$counts_cache
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
            c(
              "x" = "The Proximity score table does not contain the required {.str join_count} column"
            )
          )
        }
        if (!"join_count_expected_mean" %in% colnames(proximity_scores)) {
          cli::cli_abort(
            c(
              "x" = "The Proximity score table does not contain the required {.str join_count_expected_mean} column"
            )
          )
        }
        proximity_scores <- proximity_scores %>%
          mutate(
            log2_ratio = log2(
              pmax(join_count, 1) / pmax(join_count_expected_mean, 1)
            )
          )
      }

      proximity_scores <- proximity_scores %>% compute()

      if (!lazy) {
        proximity_scores <- proximity_scores %>% as_tibble()
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
    #' @param umi_data_type One of "int64" or "string". Default is "int64".
    #'  - "int64": The UMIs are encoded as int64
    #'  - "string": The UMIs are encoded as character
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
                                   umi_data_type = c("int64", "string"),
                                   lazy = FALSE,
                                   include_all_columns = FALSE) {
      self$check_connection()
      assert_vector(components, "character", n = 1, allow_null = TRUE)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)
      assert_single_value(lazy, type = "bool")
      assert_single_value(include_all_columns, type = "bool")

      umi_data_type <- match.arg(
        umi_data_type,
        choices = c("int64", "string")
      )

      if (!lazy && !is.null(components)) {
        el <- tbl(
          private$con,
          dbplyr::sql(private$build_component_source("edgelist", components))
        )
      } else {
        el <- tbl(private$con, "edgelist")
        if (!is.null(components)) {
          el <- el %>%
            filter(component %in% components)
        }
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

      if (lazy) {
        # Register as a temporary VIEW instead of materializing a table
        sql_query <- dbplyr::sql_render(el)
        DBI::dbExecute(
          private$con,
          paste0(
            "CREATE OR REPLACE TEMPORARY VIEW edgelist_modified AS ",
            sql_query
          )
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
    components_layout = function(components,
                                 add_marker_counts = FALSE,
                                 verbose = TRUE) {
      self$check_connection()
      assert_vector(components, "character", allow_null = TRUE, n = 1)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)
      assert_single_value(add_marker_counts, "bool")

      # Construct the SQL query
      base_select <- "SELECT CAST(index AS STRING) AS name, x, y, z FROM layouts"

      components <- components %||% (self$counts() %>% colnames())

      if (verbose) {
        cli::cli_alert_info(
          "Fetching {.val {length(components)}} component layouts..."
        )
      }

      use_index <- length(components) > 1L || !is.null(private$component_index_cache[["layouts"]])
      idx <- if (use_index) private$get_component_index("layouts") else NULL
      contiguous <- !is.null(idx) && isTRUE(attr(idx, "contiguous"))
      sub_idx <- if (!is.null(idx)) idx[match(components, idx$component), , drop = FALSE] else NULL

      lapply_func <- ifelse(verbose, pbapply::pblapply, lapply)

      layout_list <- lapply_func(seq_along(components), function(i) {
        component <- components[[i]]
        where_sql <- if (!is.null(sub_idx) && !is.na(sub_idx$min_id[[i]])) {
          if (contiguous) {
            sprintf(" WHERE rowid BETWEEN %.0f AND %.0f", sub_idx$min_id[[i]], sub_idx$max_id[[i]])
          } else {
            sprintf(
              " WHERE (rowid BETWEEN %.0f AND %.0f) AND component = '%s'",
              sub_idx$min_id[[i]], sub_idx$max_id[[i]], component
            )
          }
        } else {
          glue::glue(" WHERE component IN ('{component}')")
        }
        DBI::dbGetQuery(private$con, paste0(base_select, where_sql)) %>%
          as_tibble()
      }) %>%
        set_names(nm = components)

      if (add_marker_counts) {
        if (verbose) {
          cli::cli_alert_info("Fetching marker counts...")
        }
        marker_counts <- self$components_marker_counts(components)
        if (verbose) {
          cli::cli_alert_info("Adding marker counts to layout tables...")
        }
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
    components_marker_counts = function(components,
                                        as_sparse = FALSE,
                                        verbose = FALSE) {
      self$check_connection()
      assert_vector(components, "character", allow_null = TRUE, n = 1)
      assert_x_in_y(components, self$counts() %>% colnames(), allow_null = TRUE)

      components <- components %||% (self$counts() %>% colnames())

      use_index <- length(components) > 1L || !is.null(private$component_index_cache[["edgelist"]])
      idx <- if (use_index) private$get_component_index("edgelist") else NULL
      contiguous <- !is.null(idx) && isTRUE(attr(idx, "contiguous"))
      sub_idx <- if (!is.null(idx)) idx[match(components, idx$component), , drop = FALSE] else NULL

      lapply_func <- ifelse(verbose, pbapply::pblapply, lapply)

      component_list <- lapply_func(seq_along(components), function(i) {
        x <- components[[i]]
        where_sql <- if (!is.null(sub_idx) && !is.na(sub_idx$min_id[[i]])) {
          if (contiguous) {
            sprintf("WHERE rowid BETWEEN %.0f AND %.0f", sub_idx$min_id[[i]], sub_idx$max_id[[i]])
          } else {
            sprintf(
              "WHERE (rowid BETWEEN %.0f AND %.0f) AND component = '%s'",
              sub_idx$min_id[[i]], sub_idx$max_id[[i]], x
            )
          }
        } else {
          sprintf("WHERE component IN ('%s')", x)
        }
        sql_query <- paste0(
          "WITH src AS (SELECT umi1, umi2, marker_1, marker_2 FROM edgelist ", where_sql, ") ",
          "SELECT CAST(umi1 AS STRING) AS name, marker_1 AS marker FROM src ",
          "UNION ",
          "SELECT CAST(umi2 AS STRING) AS name, marker_2 AS marker FROM src"
        )
        data <- DBI::dbGetQuery(private$con, sql_query)
        name_f <- factor(data$name, levels = unique(data$name))
        marker_f <- factor(data$marker, levels = unique(data$marker))
        m <- Matrix::sparseMatrix(
          i = as.integer(name_f),
          j = as.integer(marker_f),
          x = 1,
        )
        rownames(m) <- levels(name_f)
        colnames(m) <- levels(marker_f)
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
    #' Load cell graphs for components as a CellGraphList
    #'
    #' This method queries the \code{edgelist} table to construct bipartite
    #' \code{CellGraph} objects for the specified \code{components}. Each graph
    #' is converted to a \code{tbl_graph} and assigned attributes \code{type = "bipartite"},
    #' \code{component_id}, and \code{assay_type = "PNA"}. Node markers are
    #' returned as a sparse count matrix.
    #'
    #' @param components A character vector of component IDs.
    #' @param markers An optional character vector of all possible marker names to set fixed
    #'   columns on the counts matrix. If \code{NULL}, marker names are derived from the data.
    #' @param add_layouts A logical specifying whether to load pre-computed 3D layouts
    #'   from the \code{layouts} table into each \code{CellGraph}. Default is \code{FALSE}.
    #' @param verbose Print progress messages.
    #'
    #' @examples
    #' # Connect to a PNA database and load cell graphs
    #' pxl_file <- minimal_pna_pxl_file()
    #' db <- PixelDB$new(pxl_file)
    #' cgl <- db$components_cell_graph(c("0a45497c6bfbfb22"))
    #' cgl
    #'
    #' @return A \code{CellGraphList} object containing \code{CellGraph} objects for the components.
    #'
    components_cell_graph = function(components,
                                     markers = NULL,
<<<<<<< HEAD
                                     add_layouts = FALSE,
=======
>>>>>>> 8a4ef9ac (Added marginal_saturation methods)
                                     verbose = FALSE) {
      self$check_connection()
      assert_vector(components, "character", n = 1)
      assert_unique(components)
      assert_x_in_y(components, colnames(self$counts()))
      assert_vector(markers, "character", allow_null = TRUE, n = 1)
      assert_single_value(add_layouts, "bool")
      assert_single_value(verbose, "bool")

      if (add_layouts && !"layouts" %in% self$names()) {
        cli::cli_abort(c("x" = "Layouts are missing from the PXL file."))
      }

      # Temporarily disable DuckDB's internal query progress bar to prevent CLI display conflict
      orig_progress_bar <- DBI::dbGetQuery(
        private$con,
        "SELECT current_setting('enable_progress_bar') AS val"
      )$val[[1]]
      DBI::dbExecute(private$con, "SET enable_progress_bar = false;")
      on.exit(
        {
          if (DBI::dbIsValid(private$con)) {
            DBI::dbExecute(
              private$con,
              glue::glue("SET enable_progress_bar = {orig_progress_bar};")
            )
          }
        },
        add = TRUE
      )

      use_el_index <- length(components) > 1L || !is.null(private$component_index_cache[["edgelist"]])
      idx_el <- if (use_el_index) private$get_component_index("edgelist") else NULL
      el_contig <- !is.null(idx_el) && isTRUE(attr(idx_el, "contiguous"))
      sub_el <- if (!is.null(idx_el)) idx_el[match(components, idx_el$component), , drop = FALSE] else NULL

      idx_lay <- NULL
      lay_contig <- FALSE
      sub_lay <- NULL
      if (add_layouts) {
        use_lay_index <- length(components) > 1L || !is.null(private$component_index_cache[["layouts"]])
        idx_lay <- if (use_lay_index) private$get_component_index("layouts") else NULL
        lay_contig <- !is.null(idx_lay) && isTRUE(attr(idx_lay, "contiguous"))
        sub_lay <- if (!is.null(idx_lay)) idx_lay[match(components, idx_lay$component), , drop = FALSE] else NULL
      }

      if (verbose) {
        pb <- cli::cli_progress_bar(
          name = "Loading CellGraphs",
          total = length(components),
          clear = FALSE
        )
      }

      cg_list <- lapply(seq_along(components), function(i) {
        cell <- components[[i]]
        el_where <- if (!is.null(sub_el) && !is.na(sub_el$min_id[[i]])) {
          if (el_contig) {
            sprintf("WHERE rowid BETWEEN %.0f AND %.0f", sub_el$min_id[[i]], sub_el$max_id[[i]])
          } else {
            sprintf(
              "WHERE (rowid BETWEEN %.0f AND %.0f) AND component = '%s'",
              sub_el$min_id[[i]], sub_el$max_id[[i]], cell
            )
          }
        } else {
          sprintf("WHERE component = '%s'", cell)
        }

        e <- DBI::dbGetQuery(
          private$con,
          paste0("SELECT umi1, umi2, marker_1, marker_2 FROM edgelist ", el_where)
        )

        E <- nrow(e)
        if (E == 0L) {
          cli::cli_abort(
            c(
              "x" = "No edges found for component {.val {cell}}"
            )
          )
        }

        all_umi <- c(e$umi1, e$umi2)
        dup <- duplicated(all_umi)
        u <- all_umi[!dup]
        n <- length(u)

        # Fast 64-bit integer indexing using bit64 hashmap
        hu <- bit64::hashmap(u)
        idx <- bit64::hashpos(hu, all_umi)

        edge_mat <- matrix(idx, ncol = 2, byrow = FALSE)
        g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)

        # Convert int64 to character for vertex names and determine node_type in one pass
        u_names <- as.character(u)
        node_type <- rep.int("umi2", n)
        node_type[which(!dup) <= E] <- "umi1"
        igraph::vertex_attr(g) <- list(name = u_names, node_type = node_type)

        # Convert to tbl_graph and set bipartite attributes matching .load_pna_as_bipartite
        tbl_g <- tidygraph::as_tbl_graph(g, directed = FALSE)
        attr(tbl_g, "type") <- "bipartite"
        attr(tbl_g, "component_id") <- cell
        attr(tbl_g, "assay_type") <- "PNA"

        # Subset markers using logical vector (!dup) to avoid matching
        marker_char <- c(e$marker_1, e$marker_2)[!dup]
        if (is.null(markers)) {
          marker_factor <- as.factor(marker_char)
          marker_levels <- levels(marker_factor)
          j <- as.integer(marker_factor)
          n_markers <- length(marker_levels)
        } else {
          j <- match(marker_char, markers)
          marker_levels <- markers
          n_markers <- length(markers)
        }

        X <- Matrix::sparseMatrix(
          i = seq_len(n),
          j = j,
          x = 1,
          dims = c(n, n_markers),
          dimnames = list(NULL, marker_levels)
        )

        layout_slot <- NULL
        if (add_layouts) {
          lay_where <- if (!is.null(sub_lay) && !is.na(sub_lay$min_id[[i]])) {
            if (lay_contig) {
              sprintf("WHERE rowid BETWEEN %.0f AND %.0f", sub_lay$min_id[[i]], sub_lay$max_id[[i]])
            } else {
              sprintf(
                "WHERE (rowid BETWEEN %.0f AND %.0f) AND component = '%s'",
                sub_lay$min_id[[i]], sub_lay$max_id[[i]], cell
              )
            }
          } else {
            sprintf("WHERE component = '%s'", cell)
          }
          lay_df <- DBI::dbGetQuery(
            private$con,
            paste0("SELECT index, x, y, z FROM layouts ", lay_where)
          )
          h_lay <- bit64::hashmap(lay_df$index)
          lay_pos <- bit64::hashpos(h_lay, u)
          layout_slot <- list(
            wpmds_3d = data.frame(
              x = lay_df$x[lay_pos],
              y = lay_df$y[lay_pos],
              z = lay_df$z[lay_pos]
            )
          )
        }

        cg <- methods::new(
          Class = "CellGraph",
          cellgraph = tbl_g,
          nodes = u_names,
          counts = X,
          layout = layout_slot,
          layers = list(),
          meta.data = .empty_node_meta(n),
          reductions = list()
        )

        if (verbose) {
          cli::cli_progress_update(id = pb)
        }

        cg
      })

      if (verbose) {
        cli::cli_progress_done(id = pb)
      }

      names(cg_list) <- components
      CreateCellGraphList(cg_list)
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
      table_name <- match.arg(
        table_name,
        choices = c("proximity", "edgelist", "layouts")
      )
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
    counts_cache = NULL,
    component_index_cache = NULL,
    get_component_index = function(table_name) {
      if (is.null(private$component_index_cache[[table_name]])) {
        idx <- DBI::dbGetQuery(
          private$con,
          paste0(
            "SELECT component, MIN(rowid) AS min_id, MAX(rowid) AS max_id, COUNT(*) AS n_rows ",
            "FROM ", table_name, " GROUP BY component"
          )
        )
        idx$component <- as.character(idx$component)
        idx$min_id <- as.numeric(idx$min_id)
        idx$max_id <- as.numeric(idx$max_id)
        idx$n_rows <- as.numeric(idx$n_rows)
        attr(idx, "contiguous") <- all((idx$max_id - idx$min_id + 1) == idx$n_rows)
        private$component_index_cache[[table_name]] <- idx
      }
      private$component_index_cache[[table_name]]
    },
    build_component_source = function(table_name, components = NULL, select_expr = "*") {
      if (is.null(components)) {
        return(paste0("SELECT ", select_expr, " FROM ", table_name))
      }
      idx <- private$get_component_index(table_name)
      if (!isTRUE(attr(idx, "contiguous"))) {
        in_sql <- glue::glue_sql("{components*}", .con = private$con)
        return(paste0(
          "SELECT ", select_expr, " FROM ", table_name,
          " WHERE component IN (", in_sql, ")"
        ))
      }
      sub_idx <- idx[idx$component %in% components, , drop = FALSE]
      if (nrow(sub_idx) == 0L) {
        return(paste0("SELECT ", select_expr, " FROM ", table_name, " WHERE 1 = 0"))
      }
      sub_idx <- sub_idx[order(sub_idx$min_id), , drop = FALSE]
      mins <- sub_idx$min_id
      maxs <- sub_idx$max_id
      n <- nrow(sub_idx)

      branches <- character()
      merged_min <- mins[[1]]
      merged_max <- maxs[[1]]
      if (n > 1L) {
        for (i in 2:n) {
          if (mins[[i]] <= merged_max + 1) {
            merged_max <- max(merged_max, maxs[[i]])
          } else {
            branches <- c(
              branches,
              sprintf(
                "SELECT %s FROM %s WHERE rowid BETWEEN %.0f AND %.0f",
                select_expr, table_name, merged_min, merged_max
              )
            )
            merged_min <- mins[[i]]
            merged_max <- maxs[[i]]
          }
        }
      }
      branches <- c(
        branches,
        sprintf(
          "SELECT %s FROM %s WHERE rowid BETWEEN %.0f AND %.0f",
          select_expr, table_name, merged_min, merged_max
        )
      )
      paste(branches, collapse = " UNION ALL ")
    },
    finalize = function() {
      if (DBI::dbIsValid(private$con)) {
        DBI::dbDisconnect(private$con)
      }
    }
  ),
  cloneable = FALSE
)
