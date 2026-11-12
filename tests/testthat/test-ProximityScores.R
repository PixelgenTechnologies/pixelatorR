library(dplyr)

pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expect_no_error(seur_obj <- suppressWarnings(ReadPNA_Seurat(pxl_file, verbose = FALSE)))
  expect_no_error(pna_assay <- seur_obj[["PNA"]])
  expect_no_error(proximity <- ReadPNA_proximity(pxl_file, verbose = FALSE, calc_log2_ratio = TRUE))

  test_that("ProximityScores getter/setter works as expected", {
    # Seurat
    expect_no_error(proximity_from_seur <- ProximityScores(seur_obj))
    expect_identical(proximity_from_seur, proximity)
    expect_no_error(proximity_from_seur <- ProximityScores(seur_obj, add_marker_counts = TRUE))
    expect_identical(select(proximity_from_seur, -matches("count_\\d")), proximity)
    expect_no_error(proximity_from_seur <- ProximityScores(seur_obj, add_marker_proportions = TRUE))
    expect_no_error(ProximityScores(seur_obj) <- NULL)
    expect_length(ProximityScores(seur_obj), 0)

    # PNAAssay
    expect_no_error(ProximityScores(pna_assay, add_marker_counts = TRUE))
    expect_no_error(ProximityScores(pna_assay) <- NULL)
    expect_length(ProximityScores(pna_assay), 0)

    # Lazy load
    expect_no_error(proximity_from_seur_lazy <- ProximityScores(seur_obj, lazy = TRUE, add_marker_counts = TRUE))
    expect_no_error(proximity_from_seur_lazy <- ProximityScores(seur_obj, lazy = TRUE, add_marker_proportions = TRUE))
    expect_true(all(c("count_1", "count_2", "p1", "p2") %in% colnames(proximity_from_seur_lazy)))
    expect_identical(nrow(collect(proximity_from_seur_lazy)), nrow(proximity))
    expect_true(all(colnames(proximity) %in% colnames(proximity_from_seur_lazy)))

    DBI::dbDisconnect(proximity_from_seur_lazy$src$con)
  })

  test_that("ProximityScores fails with invalid input", {
    expect_error(ProximityScores("Invalid"))
    expect_error(ProximityScores(seur_obj, assay = FALSE))
    expect_error(ProximityScores(seur_obj, meta_data_columns = "Invalid"))
    expect_error(ProximityScores(seur_obj, add_marker_counts = "Invalid"))
    expect_error(ProximityScores(seur_obj, add_marker_proportions = "Invalid"))
    expect_error(ProximityScores(seur_obj) <- "Invalid")
  })
}

test_that("join_proximity_metadata.tbl_lazy handles join_by, NULL by, and duplicate keys", {
  seur_lazy <- suppressWarnings(ReadPNA_Seurat(pxl_file, load_proximity_scores = FALSE, verbose = FALSE))
  prox <- ProximityScores(seur_lazy, lazy = TRUE)
  on.exit(DBI::dbDisconnect(prox$src$con, shutdown = TRUE), add = TRUE)
  expect_true(isTRUE(attr(prox, "pxl_unique_keys")))

  # join_by with renamed key preserves pxl_unique_keys when unique
  m_uniq <- tibble::tibble(cell_id = colnames(seur_lazy), grp = seq_len(ncol(seur_lazy)))
  j_jb <- join_proximity_metadata(prox, m_uniq, by = dplyr::join_by(component == cell_id))
  expect_true(isTRUE(attr(j_jb, "pxl_unique_keys")))

  # by = NULL natural join preserves pxl_unique_keys when unique
  m_nat <- tibble::tibble(component = colnames(seur_lazy), grp = seq_len(ncol(seur_lazy)))
  j_null <- suppressMessages(join_proximity_metadata(prox, m_nat, by = NULL))
  expect_true(isTRUE(attr(j_null, "pxl_unique_keys")))

  # Duplicate RHS keys drop pxl_unique_keys
  m_dup <- tibble::tibble(cell_id = rep(colnames(seur_lazy)[1], 2), grp = c("a", "b"))
  j_dup <- join_proximity_metadata(prox, m_dup, by = dplyr::join_by(component == cell_id))
  expect_null(attr(j_dup, "pxl_unique_keys"))
})

