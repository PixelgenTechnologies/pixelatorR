pxl_file <- minimal_pna_pxl_file()

test_that("PixelDB initialize methods works as expected", {
  expect_no_error(db <- PixelDB$new(pxl_file))
  expect_true(inherits(db, "R6"))
  con <- db$.__enclos_env__$private$con
  expect_true(DBI::dbIsValid(con))
  db$close()
  expect_false(DBI::dbIsValid(con))
})

test_that("PixelDB finalizer eventually closes connections", {
  db <- PixelDB$new(pxl_file)
  con <- db$.__enclos_env__$private$con
  rm(db)
  gc()

  # Instead of strict expectation:
  tryCatch(
    expect_false(DBI::dbIsValid(con)),
    error = function(e) skip("Finalizer timing is platform-dependent")
  )
})

test_that("PixelDB computes proximity scores in a temporary table", {
  db <- PixelDB$new(pxl_file)
  on.exit(db$close())

  proximity <- db$compute_proximity_scores(
    components = "0a45497c6bfbfb22",
    markers = c("B2M", "HLA-ABC"),
    name = "computed_proximity"
  )

  expect_s3_class(proximity, "tbl_lazy")
  expect_setequal(
    colnames(proximity),
    c(
      "component",
      "marker_1",
      "marker_2",
      "join_count",
      "join_count_expected_mean",
      "log2_ratio"
    )
  )

  result <- proximity %>% collect()
  expect_gt(nrow(result), 0)
  expect_true(all(result$component == "0a45497c6bfbfb22"))
  expect_true(all(result$marker_1 %in% c("B2M", "HLA-ABC")))
  expect_true(all(result$marker_2 %in% c("B2M", "HLA-ABC")))
  expect_true(all(result$marker_1 <= result$marker_2))
  expect_true(all(is.finite(result$log2_ratio)))

  proximity_z <- db$compute_proximity_scores(
    components = "0a45497c6bfbfb22",
    markers = c("B2M", "HLA-ABC"),
    calc_z_score = TRUE,
    name = "computed_proximity_z"
  )
  expect_true(all(c("join_count_expected_sd", "join_count_z", "join_count_p") %in% colnames(proximity_z)))
  result_z <- proximity_z %>% collect()
  expect_true(all(is.finite(result_z$join_count_z)))

  # Verify binomial SD formula sqrt(n * p * (1 - p)) for both homotypic and heterotypic pairs
  n_edges <- db$components_edgelist("0a45497c6bfbfb22") %>% nrow()
  p <- result_z$join_count_expected_mean / n_edges
  expect_equal(result_z$join_count_expected_sd, sqrt(n_edges * p * (1 - p)))

  table_info <- db$info() %>%
    filter(name == "computed_proximity")
  expect_equal(nrow(table_info), 1)
  expect_true(table_info$temporary)
})

test_that("PixelDB methods work as expected", {
  # close method
  expect_no_error(db <- PixelDB$new(pxl_file))
  con <- db$.__enclos_env__$private$con
  expect_no_error(db$close())
  if (.Platform$OS.type == "unix") {
    expect_true(!DBI::dbIsValid(con))
  }

  # info method
  expect_no_error(db <- PixelDB$new(pxl_file))
  expect_s3_class(db$info(), "tbl_df")
  expect_equal(
    names(db$info()),
    c(
      "database",
      "schema",
      "name",
      "column_names",
      "column_types",
      "temporary"
    )
  )
  expect_equal(dim(db$info()), c(10, 6))

  # reconnect method
  con <- db$.__enclos_env__$private$con
  DBI::dbDisconnect(con)
  expect_no_error(db$reconnect())

  # query method
  expect_no_error(mdata <- db$query("SELECT * FROM metadata"))
  expected_mdata <-
    structure(
      list(value = "{\"sample_name\":\"PNA055_Sample07_S7\",\"version\":\"0.1.0\",\"technology\":\"single-cell-pna\",\"panel_name\":\"pna-rnd-158plex-final\",\"panel_version\":\"0.1.0\",\"post_analysis\":{\"params\":{\"degree-analysis-per-marker-type\":{\"degree-analysis-per-marker-type\":{}}}}}"),
      class = "data.frame",
      row.names = c(NA, -1L)
    )
  expect_equal(mdata, expected_mdata)

  # names method
  expect_no_error(table_names <- db$names())
  expect_equal(
    table_names,
    c(
      "__adata__X",
      "__adata__obs",
      "__adata__obsm_clr",
      "__adata__obsm_log1p",
      "__adata__uns",
      "__adata__var",
      "edgelist",
      "layouts",
      "metadata",
      "proximity"
    )
  )

  # fetch_table method
  expect_no_error(X <- db$fetch_table("__adata__X"))
  expect_equal(dim(X), c(5, 159))
  expect_s3_class(X, "data.frame")

  # fetch_table_subset method
  expect_no_error(prox <- db$fetch_table_subset("proximity", columns_filter = list("component" = "0a45497c6bfbfb22", "marker_1" = "B2M", "marker_2" = "B2M")))
  expect_equal(dim(prox), c(1, 8))
  expect_s3_class(prox, "data.frame")

  # counts method
  expect_no_error(X <- db$counts())
  expected_X <-
    new(
      "dgCMatrix",
      i = c(0L, 1L, 0L, 1L),
      p = c(0L, 2L, 4L),
      Dim = c(
        2L,
        2L
      ),
      Dimnames = list(
        c("HLA-ABC", "B2M"),
        c(
          "0a45497c6bfbfb22",
          "2708240b908e2eba"
        )
      ),
      x = c(865, 1182, 2077, 3448),
      factors = list()
    )
  expect_equal(X[1:2, 1:2], expected_X)

  # proximity method
  expect_no_error(prox <- db$proximity())
  expected_prox <-
    structure(
      list(
        marker_1 = c("CD56", "CD56"),
        marker_2 = c(
          "CD56",
          "mIgG2b"
        ),
        join_count = c(0, 0),
        join_count_expected_mean = c(
          0,
          0.03
        ),
        join_count_expected_sd = c(0, 0.171446607997765),
        join_count_z = c(0, -0.03),
        join_count_p = c(0.5, 0.488033526585887),
        component = c(
          "c3c393e9a17c1981",
          "c3c393e9a17c1981"
        ),
        log2_ratio = c(0, 0)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  expect_equal(prox |> head(n = 2), expected_prox)

  # cell_meta method
  expect_no_error(cell_meta <- db$cell_meta())
  expect_equal(dim(cell_meta), c(5, 23))
  expect_s3_class(cell_meta, "data.frame")

  # protein_meta method
  expect_no_error(protein_meta <- db$protein_meta())
  expect_equal(dim(protein_meta), c(158, 5))
  expect_s3_class(protein_meta, "data.frame")

  # run_meta method
  expect_no_error(run_meta <- db$run_meta())
  expect_equal(dim(run_meta), c(1, 6))
  expect_s3_class(run_meta, "data.frame")

  # components_edgelist method
  expect_no_error(el <- db$components_edgelist("c3c393e9a17c1981"))
  expect_equal(
    el |> sapply(class),
    c(
      marker_1 = "character",
      marker_2 = "character",
      umi1 = "integer64",
      umi2 = "integer64",
      component = "character"
    )
  )
  expect_equal(dim(el), c(110657, 5))

  expect_no_error(el <- db$components_edgelist("c3c393e9a17c1981", umi_data_type = "string"))
  expect_equal(
    el |> sapply(class),
    c(
      marker_1 = "character",
      marker_2 = "character",
      umi1 = "character",
      umi2 = "character",
      component = "character"
    )
  )
  expect_equal(dim(el), c(110657, 5))

  expect_no_error(el_lazy <- db$components_edgelist(c("0a45497c6bfbfb22", "c3c393e9a17c1981"), lazy = TRUE))
  expect_no_error(el <- db$components_edgelist(c("0a45497c6bfbfb22", "c3c393e9a17c1981"), lazy = FALSE))

  expect_s3_class(
    el_lazy, "tbl_lazy"
  )

  expect_equal(
    el_lazy |> collect(),
    el
  )

  # components_marker_counts method
  expect_no_error(mc <- db$components_marker_counts("0a45497c6bfbfb22"))
  expect_type(mc, "list")
  expect_equal(names(mc), "0a45497c6bfbfb22")
  expect_equal(dim(mc[[1]]), c(43543, 150))

  # components_cell_graph method (single and multiple components, markers, and add_layouts)
  cells_two <- c("0a45497c6bfbfb22", "2708240b908e2eba")
  expect_no_error(cgl_one <- db$components_cell_graph(cells_two[1]))
  expect_s3_class(cgl_one, "CellGraphList")
  expect_named(cgl_one, cells_two[1])
  expect_s4_class(cgl_one[[1]], "CellGraph")
  expect_equal(attr(cgl_one[[1]]@cellgraph, "type"), "bipartite")
  expect_equal(attr(cgl_one[[1]]@cellgraph, "component_id"), cells_two[1])
  expect_equal(attr(cgl_one[[1]]@cellgraph, "assay_type"), "PNA")
  expect_null(cgl_one[[1]]@layout)

  # Multiple components trigger the cached rowid component index
  expect_no_error(
    cgl_two <- db$components_cell_graph(
      cells_two,
      markers = rownames(db$counts()),
      add_layouts = TRUE
    )
  )
  expect_s3_class(cgl_two, "CellGraphList")
  expect_named(cgl_two, cells_two)
  expect_equal(cgl_two[[1]]@nodes, cgl_one[[1]]@nodes)
  expect_equal(colnames(cgl_two[[1]]@counts), rownames(db$counts()))
  expect_named(cgl_two[[1]]@layout, "wpmds_3d")
  expect_equal(nrow(cgl_two[[1]]@layout$wpmds_3d), length(cgl_two[[1]]@nodes))
  expect_equal(colnames(cgl_two[[1]]@layout$wpmds_3d), c("x", "y", "z"))

  # Subset of markers (and duplicate markers) preserves node count and filters non-matching nodes to 0
  subset_markers <- c("HLA-ABC", "B2M", "HLA-ABC")
  expect_no_error(
    cgl_sub <- db$components_cell_graph(cells_two[1], markers = subset_markers)
  )
  expect_equal(
    dim(cgl_sub[[1]]@counts),
    c(length(cgl_one[[1]]@nodes), 2L)
  )
  expect_equal(colnames(cgl_sub[[1]]@counts), c("HLA-ABC", "B2M"))
  expect_equal(
    cgl_sub[[1]]@counts,
    cgl_one[[1]]@counts[, c("HLA-ABC", "B2M"), drop = FALSE]
  )

  # Sorted edgelist skips the full-table rowid index scan, while unsorted layouts falls back to rowid index
  sorted_cache <- db$.__enclos_env__$private$table_sorted_cache
  expect_true(isTRUE(sorted_cache$edgelist))
  expect_false(isTRUE(sorted_cache$layouts))
  idx_cache <- db$.__enclos_env__$private$component_index_cache
  expect_false("edgelist" %in% names(idx_cache))
  expect_true("layouts" %in% names(idx_cache))
  expect_equal(nrow(idx_cache$layouts), ncol(db$counts()))

  # export parquet method
  tmp_parquet_file <- fs::file_temp(ext = "parquet")
  expect_no_error(db$export_parquet(tmp_parquet_file, table_name = "proximity"))
  expect_true(fs::file_exists(tmp_parquet_file))
  tmp_parquet_file <- fs::file_temp(ext = "parquet")
  expect_no_error(db$export_parquet(tmp_parquet_file, table_name = "edgelist"))
  expect_true(fs::file_exists(tmp_parquet_file))
  tmp_parquet_file <- fs::file_temp(ext = "parquet")
  expect_no_error(db$export_parquet(tmp_parquet_file, table_name = "layouts"))
  expect_true(fs::file_exists(tmp_parquet_file))

  expect_no_error(db$close())
})


test_that("PixelDB methods fails with invalid input", {
  expect_error(suppressWarnings(db <- PixelDB$new("Invalid")))
  expect_no_error(db <- PixelDB$new(pxl_file))
  expect_error(db$query("Invalid"))
  expect_error(db$fetch_table("Invalid"))
  expect_error(db$fetch_table_subset("Invalid"))
  expect_error(db$proximity(calc_log2_ratio = "Invalid"))
  expect_error(db$compute_proximity_scores(components = 1))
  expect_error(db$compute_proximity_scores(markers = 1))
  expect_error(db$compute_proximity_scores(name = c("one", "two")))
  expect_error(db$compute_proximity_scores(calc_z_score = "Invalid"))
  expect_error(db$compute_proximity_scores(min_marker_count = -1L))
  expect_error(db$compute_proximity_scores(batch_size = 0L))
  expect_error(db$components_edgelist("Invalid"))
  expect_error(db$components_edgelist("3898b03349c6e28d", umi_data_type = "Invalid"))
  expect_error(db$components_edgelist("3898b03349c6e28d", include_all_columns = "Invalid"))
  expect_error(db$components_layout("Invalid"))
  expect_error(db$components_layout("3898b03349c6e28d", add_marker_counts = "Invalid"))
  expect_error(db$components_layout("3898b03349c6e28d", verbose = FALSE))
  expect_error(db$components_marker_counts("Invalid"))
  expect_error(db$components_cell_graph("Invalid"))
  expect_error(db$components_cell_graph(c("0a45497c6bfbfb22", "0a45497c6bfbfb22")))
  expect_error(db$components_cell_graph("0a45497c6bfbfb22", add_layouts = "Invalid"))
  expect_no_error(db$close())
})

test_that("PixelDB components_cell_graph falls back to rowid index when edgelist is unsorted", {
  tmp_pxl <- tempfile(fileext = ".pxl")
  fs::file_copy(pxl_file, tmp_pxl)
  on.exit(unlink(tmp_pxl), add = TRUE)

  # Reverse component order in edgelist so it is no longer sorted ascending
  con_rw <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp_pxl, read_only = FALSE)
  DBI::dbExecute(con_rw, "CREATE OR REPLACE TABLE edgelist AS SELECT * FROM edgelist ORDER BY component DESC")
  DBI::dbDisconnect(con_rw, shutdown = TRUE)

  db_unsorted <- PixelDB$new(tmp_pxl, threads = 1L)
  on.exit(db_unsorted$close(), add = TRUE)
  cells_two <- c("0a45497c6bfbfb22", "2708240b908e2eba")
  cgl <- db_unsorted$components_cell_graph(cells_two)
  expect_named(cgl, cells_two)
  expect_false(isTRUE(db_unsorted$.__enclos_env__$private$table_sorted_cache$edgelist))
  expect_true("edgelist" %in% names(db_unsorted$.__enclos_env__$private$component_index_cache))
})

test_that("PixelDB respects pixelatorR.read_only option", {
  tmp_pxl <- tempfile(fileext = ".pxl")
  fs::file_copy(pxl_file, tmp_pxl)
  on.exit(unlink(tmp_pxl), add = TRUE)

  # Default is read-only
  db_ro <- PixelDB$new(tmp_pxl)
  on.exit(db_ro$close(), add = TRUE)
  expect_true(db_ro$.__enclos_env__$private$read_only)
  expect_error(
    DBI::dbExecute(db_ro$.__enclos_env__$private$con, "CREATE TABLE tmp_test (x INT)")
  )
  db_ro$close()

  # Invalid option value fails validation
  old_opt <- options(pixelatorR.read_only = "invalid")
  on.exit(options(old_opt), add = TRUE)
  expect_error(PixelDB$new(tmp_pxl))
  options(old_opt)

  # Setting pixelatorR.read_only = FALSE opens writable connection and persists across reconnect
  old_opt <- options(pixelatorR.read_only = FALSE)
  db_rw <- PixelDB$new(tmp_pxl)
  options(old_opt)
  on.exit(db_rw$close(), add = TRUE)

  expect_false(db_rw$.__enclos_env__$private$read_only)
  expect_no_error(
    DBI::dbExecute(db_rw$.__enclos_env__$private$con, "CREATE TABLE tmp_test (x INT)")
  )
  expect_message(db_rw$reconnect(), "Connected")
  expect_no_error(
    DBI::dbExecute(db_rw$.__enclos_env__$private$con, "DROP TABLE tmp_test")
  )
})

write_null_pxl <- function(path, reason) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # VARCHAR stores the metadata JSON as text. A JSON column needs the DuckDB
  # json extension, which R CMD check cannot autoinstall.
  DBI::dbExecute(con, "CREATE TABLE metadata (value VARCHAR)")
  payload <- list(
    sample_name = "null-sample",
    version = "0.30.0",
    null = TRUE
  )
  if (!is.null(reason)) {
    payload$null_reason <- reason
  }
  meta_json <- as.character(jsonlite::toJSON(payload, auto_unbox = TRUE))
  DBI::dbExecute(
    con,
    "INSERT INTO metadata (value) VALUES (?)",
    params = list(meta_json)
  )
}

test_that("PixelDB throws descriptive error when PXL file is empty (marked as null)", {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  tmp_pxl <- file.path(tmp_dir, "null_sample.pxl")
  write_null_pxl(tmp_pxl, "Processing failed at step X")

  expect_error(
    PixelDB$new(tmp_pxl),
    class = "pixeldb_null_error",
    regexp = "The PXL file is empty \\(marked as null\\)"
  )
  expect_error(
    PixelDB$new(tmp_pxl),
    class = "pixeldb_null_error",
    regexp = "Reason: Processing failed at step X"
  )

  # Braces in the reason must stay literal and still use pixeldb_null_error.
  tmp_pxl_braces <- file.path(tmp_dir, "null_sample_braces.pxl")
  write_null_pxl(tmp_pxl_braces, "failed {step} with } brace")
  err <- expect_error(
    PixelDB$new(tmp_pxl_braces),
    class = "pixeldb_null_error"
  )
  expect_match(
    conditionMessage(err),
    "Reason: failed {step} with } brace",
    fixed = TRUE
  )
})

write_protein_meta_test_pxl <- function(
  path,
  panels_df = NULL,
  var_has_control = FALSE
) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = path)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  var_df <- data.frame(
    index = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG", "unlisted"),
    antibody_count = bit64::as.integer64(c(30L, 1L, 2L, 2L, 1L, 7L)),
    antibody_pct = c(0.6976744, 0.0232558, 0.0465116, 0.0465116, 0.0232558, 0.1627907),
    components = bit64::as.integer64(c(2L, 1L, 1L, 2L, 1L, 2L)),
    check.names = FALSE
  )
  if (var_has_control) {
    var_df$control <- c("no", "yes", "yes", "no", "yes", "no")
  }
  DBI::dbWriteTable(con, "__adata__var", var_df)

  if (!is.null(panels_df)) {
    DBI::dbWriteTable(con, "panels", panels_df)
  }
}

test_that("PixelDB$protein_meta copies missing panel columns as expected", {
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  panels_df <- data.frame(
    row_nr = bit64::as.integer64(c(4L, 3L, 2L, 1L, 0L, 5L)),
    marker_id = c("rabbit-IgG", "mIgG2b", "mIgG2a", "mIgG1", "CD3E", "panel_only"),
    control = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    full_name = c("Rabbit IgG", "Mouse IgG2b", "Mouse IgG2a", "Mouse IgG1", "CD3 epsilon", "Panel Only"),
    alt_id = c("", "", "", "", "T3E", ""),
    uniprot_id = c("", "", "", "", "P07766", ""),
    sequence_1 = c("AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC"),
    sequence_2 = c("TTTT", "TTTG", "TTTC", "TTTA", "TTGT", "TTGG"),
    sample_hashing = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    source_id = bit64::as.integer64(c(0L, 0L, 0L, 0L, 0L, 0L)),
    check.names = FALSE
  )

  # File with panels table fills missing columns matched by marker_id
  pxl_with_panels <- file.path(tmp_dir, "with_panels.pxl")
  write_protein_meta_test_pxl(pxl_with_panels, panels_df = panels_df)
  db <- PixelDB$new(pxl_with_panels)
  pm_with_panels <- db$protein_meta()
  db$close()

  expect_equal(
    pm_with_panels,
    structure(
      list(
        antibody_count = c(30L, 1L, 2L, 2L, 1L, 7L),
        antibody_pct = c(0.6976744, 0.0232558, 0.0465116, 0.0465116, 0.0232558, 0.1627907),
        components = c(2L, 1L, 1L, 2L, 1L, 2L),
        row_nr = c(0L, 1L, 2L, 3L, 4L, NA),
        control = c(FALSE, TRUE, TRUE, TRUE, TRUE, NA),
        full_name = c("CD3 epsilon", "Mouse IgG1", "Mouse IgG2a", "Mouse IgG2b", "Rabbit IgG", NA),
        alt_id = c("T3E", "", "", "", "", NA),
        uniprot_id = c("P07766", "", "", "", "", NA),
        sequence_1 = c("AACA", "AAAT", "AAAG", "AAAC", "AAAA", NA),
        sequence_2 = c("TTGT", "TTTA", "TTTC", "TTTG", "TTTT", NA),
        sample_hashing = c(FALSE, FALSE, FALSE, FALSE, FALSE, NA),
        source_id = c(0L, 0L, 0L, 0L, 0L, NA)
      ),
      row.names = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG", "unlisted"),
      class = "data.frame"
    )
  )

  # File without panels table returns only __adata__var columns
  pxl_no_panels <- file.path(tmp_dir, "no_panels.pxl")
  write_protein_meta_test_pxl(pxl_no_panels, panels_df = NULL)
  db <- PixelDB$new(pxl_no_panels)
  pm_no_panels <- db$protein_meta()
  db$close()

  expect_equal(
    pm_no_panels,
    structure(
      list(
        antibody_count = c(30L, 1L, 2L, 2L, 1L, 7L),
        antibody_pct = c(0.6976744, 0.0232558, 0.0465116, 0.0465116, 0.0232558, 0.1627907),
        components = c(2L, 1L, 1L, 2L, 1L, 2L)
      ),
      row.names = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG", "unlisted"),
      class = "data.frame"
    )
  )

  # File that already has control on __adata__var keeps those values while adding other missing columns
  pxl_var_ctrl <- file.path(tmp_dir, "var_has_control.pxl")
  write_protein_meta_test_pxl(pxl_var_ctrl, panels_df = panels_df, var_has_control = TRUE)
  db <- PixelDB$new(pxl_var_ctrl)
  pm_var_ctrl <- db$protein_meta()
  db$close()

  expect_equal(
    pm_var_ctrl,
    structure(
      list(
        antibody_count = c(30L, 1L, 2L, 2L, 1L, 7L),
        antibody_pct = c(0.6976744, 0.0232558, 0.0465116, 0.0465116, 0.0232558, 0.1627907),
        components = c(2L, 1L, 1L, 2L, 1L, 2L),
        control = c("no", "yes", "yes", "no", "yes", "no"),
        row_nr = c(0L, 1L, 2L, 3L, 4L, NA),
        full_name = c("CD3 epsilon", "Mouse IgG1", "Mouse IgG2a", "Mouse IgG2b", "Rabbit IgG", NA),
        alt_id = c("T3E", "", "", "", "", NA),
        uniprot_id = c("P07766", "", "", "", "", NA),
        sequence_1 = c("AACA", "AAAT", "AAAG", "AAAC", "AAAA", NA),
        sequence_2 = c("TTGT", "TTTA", "TTTC", "TTTG", "TTTT", NA),
        sample_hashing = c(FALSE, FALSE, FALSE, FALSE, FALSE, NA),
        source_id = c(0L, 0L, 0L, 0L, 0L, NA)
      ),
      row.names = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG", "unlisted"),
      class = "data.frame"
    )
  )

  # File where panels has no missing columns returns __adata__var unchanged
  pxl_no_missing <- file.path(tmp_dir, "no_missing_cols.pxl")
  write_protein_meta_test_pxl(
    pxl_no_missing,
    panels_df = panels_df[, c("marker_id", "control")],
    var_has_control = TRUE
  )
  db <- PixelDB$new(pxl_no_missing)
  pm_no_missing <- db$protein_meta()
  db$close()

  expect_equal(
    pm_no_missing,
    structure(
      list(
        antibody_count = c(30L, 1L, 2L, 2L, 1L, 7L),
        antibody_pct = c(0.6976744, 0.0232558, 0.0465116, 0.0465116, 0.0232558, 0.1627907),
        components = c(2L, 1L, 1L, 2L, 1L, 2L),
        control = c("no", "yes", "yes", "no", "yes", "no")
      ),
      row.names = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG", "unlisted"),
      class = "data.frame"
    )
  )

  # Malformed panels table without marker_id errors clearly
  pxl_bad_panels <- file.path(tmp_dir, "bad_panels.pxl")
  write_protein_meta_test_pxl(
    pxl_bad_panels,
    panels_df = data.frame(control = c(TRUE, FALSE))
  )
  db <- PixelDB$new(pxl_bad_panels)
  on.exit(db$close(), add = TRUE)
  expect_error(db$protein_meta(), regexp = "marker_id")
})
