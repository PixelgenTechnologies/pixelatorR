pxl_file <- minimal_pna_pxl_file()

test_that("PixelDB initialize/finalize methods works as expected", {
  expect_no_error(db <- PixelDB$new(pxl_file))
  expect_true(inherits(db, "R6"))
  con <- db$.__enclos_env__$private$con
  expect_true(DBI::dbIsValid(con))
  expect_no_error(rm(db))
  gc(full = FALSE)
  # Connection should be garbage collected
  expect_true(!DBI::dbIsValid(con))
})

test_that("PixelDB methods work as expected", {
  # close method
  expect_no_error(db <- PixelDB$new(pxl_file))
  con <- db$.__enclos_env__$private$con
  expect_no_error(db$close())
  expect_true(!DBI::dbIsValid(con))

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
      row.names = c(NA,-1L)
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
      Dim = c(2L,
              2L),
      Dimnames = list(
        c("HLA-ABC", "B2M"),
        c("0a45497c6bfbfb22",
          "2708240b908e2eba")
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
        marker_2 = c("CD56",
                     "mIgG2b"),
        join_count = c(0, 0),
        join_count_expected_mean = c(0,
                                     0.03),
        join_count_expected_sd = c(0, 0.171446607997765),
        join_count_z = c(0,-0.03),
        join_count_p = c(0.5, 0.488033526585887),
        component = c("c3c393e9a17c1981",
                      "c3c393e9a17c1981"),
        log2_ratio = c(0, 0)
      ),
      row.names = c(NA,-2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  expect_equal(prox %>% head(n = 2), expected_prox)

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
    el %>% sapply(class),
    c(
      umi1 = "integer64",
      umi2 = "integer64",
      marker_1 = "character",
      marker_2 = "character",
      component = "character"
    )
  )
  expect_equal(dim(el), c(110657, 5))

  expect_no_error(el <- db$components_edgelist("c3c393e9a17c1981", umi_data_type = "string"))
  expect_equal(
    el %>% sapply(class),
    c(
      umi1 = "character",
      umi2 = "character",
      marker_1 = "character",
      marker_2 = "character",
      component = "character"
    )
  )
  expect_equal(dim(el), c(110657, 5))

  expect_no_error(el <- db$components_edgelist("0a45497c6bfbfb22", umi_data_type = "suffixed_string"))
  expect_equal(
    el %>% sapply(class),
    c(
      umi1 = "character",
      umi2 = "character",
      marker_1 = "character",
      marker_2 = "character",
      component = "character"
    )
  )
  expect_equal(dim(el), c(97014, 5))

  # components_marker_counts method
  expect_no_error(mc <- db$components_marker_counts("0a45497c6bfbfb22"))
  expect_type(mc, "list")
  expect_equal(names(mc), "0a45497c6bfbfb22")
  expect_equal(dim(mc[[1]]), c(43543, 150))

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
  expect_error(db$components_edgelist("Invalid"))
  expect_error(db$components_edgelist("3898b03349c6e28d", umi_data_type = "Invalid"))
  expect_error(db$components_edgelist("3898b03349c6e28d", include_all_columns = "Invalid"))
  expect_error(db$components_layout("Invalid"))
  expect_error(db$components_layout("3898b03349c6e28d", add_marker_counts = "Invalid"))
  expect_error(db$components_layout("3898b03349c6e28d", verbose = FALSE))
  expect_error(db$components_marker_counts("Invalid"))
  expect_no_error(db$close())
})
