pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  test_that("ReadPNA_Seurat works as expected", {
    # default
    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, verbose = FALSE)))
    expect_s4_class(seur_obj, "Seurat")
    expect_equal(dim(seur_obj), c(158, 5))
    expect_s4_class(seur_obj[["PNA"]], ifelse(assay_version == "v3", "PNAAssay", "PNAAssay5"))

    # QC columns plus every k_core_* column; other obs columns stay out
    db <- PixelDB$new(pxl_file)
    obs_cols <- colnames(db$cell_meta())
    db$close()
    k_core_cols <- grep("^k_core_", obs_cols, value = TRUE)
    expect_gte(length(k_core_cols), 1)
    expect_equal(
      colnames(seur_obj[[]]),
      unique(c(intersect(pixelatorR:::CELL_META_COLS, obs_cols), k_core_cols))
    )

    # no proximity scores
    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, load_proximity_scores = FALSE, verbose = FALSE)))
    expect_true(length(ProximityScores(seur_obj)) == 0)

    # new assay name
    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, assay = "myassay", verbose = FALSE)))
    expect_true("myassay" %in% SeuratObject::Assays(seur_obj))

    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, return_pna_assay = FALSE, verbose = FALSE)))
    expect_s4_class(seur_obj[["PNA"]], ifelse(assay_version == "v3", "Assay", "Assay5"))

    # Including detailed meta data
    expect_no_error(suppressWarnings(seur_obj_detailed_meta <- ReadPNA_Seurat(pxl_file, detailed_meta_data = TRUE, verbose = FALSE)))
    expect_true(ncol(seur_obj_detailed_meta[[]]) > ncol(seur_obj[[]]))
  })

  test_that("ReadPNA_Seurat fails with invalid input", {
    expect_error(ReadPNA_Seurat("Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, verbose = "Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, return_pna_assay = "Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, load_proximity_scores = "Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, assay = FALSE))
  })

  test_that("ReadPNA_Seurat fails when X collapses to a vector (1 cell)", {
    # Subset to one cell and drop to a vector so ncol(X) is NULL
    trace(ReadPNA_Seurat, tracer = quote(X <- X[, 1]), at = 14, print = FALSE)
    on.exit(untrace(ReadPNA_Seurat), add = TRUE)

    expect_error(
      ReadPNA_Seurat(pxl_file, verbose = FALSE),
      regexp = "least 2|contains only one"
    )
  })

  test_that("ReadPNA_Seurat fails when X is a 1-column matrix", {
    # Inject an error by subsetting X to have only one cell, keeping it as a matrix
    trace(ReadPNA_Seurat, tracer = quote(X <- X[, 1, drop = FALSE]), at = 14, print = FALSE)
    on.exit(untrace(ReadPNA_Seurat), add = TRUE)

    expect_error(
      ReadPNA_Seurat(pxl_file, verbose = FALSE),
      regexp = "least 2|contains only one"
    )
  })

  test_that("ReadPNA_Seurat copies missing panel columns onto feature metadata", {
    tmp_pxl <- tempfile(fileext = ".pxl")
    on.exit(unlink(tmp_pxl), add = TRUE)

    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = tmp_pxl)
    DBI::dbWriteTable(
      con,
      "__adata__X",
      data.frame(
        index = c("cell_1", "cell_2"),
        CD3E = c(10L, 20L),
        mIgG1 = c(1L, 0L),
        mIgG2a = c(0L, 2L),
        mIgG2b = c(1L, 1L),
        `rabbit-IgG` = c(0L, 1L),
        check.names = FALSE
      )
    )
    DBI::dbWriteTable(
      con,
      "__adata__obs",
      data.frame(
        index = c("cell_1", "cell_2"),
        n_umi = c(12L, 24L),
        k_core_1 = c(5L, 6L),
        check.names = FALSE
      )
    )
    DBI::dbWriteTable(
      con,
      "__adata__var",
      data.frame(
        index = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG"),
        antibody_count = bit64::as.integer64(c(30L, 1L, 2L, 2L, 1L)),
        antibody_pct = c(0.8333333, 0.0277778, 0.0555556, 0.0555556, 0.0277778),
        components = bit64::as.integer64(c(2L, 1L, 1L, 2L, 1L)),
        check.names = FALSE
      )
    )
    DBI::dbWriteTable(
      con,
      "panels",
      data.frame(
        marker_id = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG"),
        control = c(FALSE, TRUE, TRUE, TRUE, TRUE),
        full_name = c("CD3 epsilon", "Mouse IgG1", "Mouse IgG2a", "Mouse IgG2b", "Rabbit IgG"),
        check.names = FALSE
      )
    )
    DBI::dbDisconnect(con, shutdown = TRUE)

    seur_obj <- suppressWarnings(
      ReadPNA_Seurat(tmp_pxl, load_proximity_scores = FALSE, verbose = FALSE)
    )
    feature_meta <- if (assay_version == "v5") {
      seur_obj[["PNA"]][[]]
    } else {
      seur_obj[["PNA"]]@meta.features
    }

    expect_equal(
      feature_meta |>
        dplyr::select(-dplyr::any_of(c("marker", "var.features", "var.features.rank"))),
      structure(
        list(
          antibody_count = c(30L, 1L, 2L, 2L, 1L),
          antibody_pct = c(0.8333333, 0.0277778, 0.0555556, 0.0555556, 0.0277778),
          components = c(2L, 1L, 1L, 2L, 1L),
          control = c(FALSE, TRUE, TRUE, TRUE, TRUE),
          full_name = c("CD3 epsilon", "Mouse IgG1", "Mouse IgG2a", "Mouse IgG2b", "Rabbit IgG")
        ),
        row.names = c("CD3E", "mIgG1", "mIgG2a", "mIgG2b", "rabbit-IgG"),
        class = "data.frame"
      )
    )
  })
}
