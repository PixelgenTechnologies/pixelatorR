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
    expect_no_error(ProximityScores(seur_obj) <- NULL)
    expect_length(ProximityScores(seur_obj), 0)

    # PNAAssay
    expect_no_error(ProximityScores(pna_assay, add_marker_counts = TRUE))
    expect_no_error(ProximityScores(pna_assay) <- NULL)
    expect_length(ProximityScores(pna_assay), 0)

    # Lazy load
    expect_no_error(proximity_from_seur_lazy <- ProximityScores(seur_obj, lazy = TRUE, add_marker_counts = TRUE))
    expect_true(all(c("count_1", "count_2") %in% colnames(proximity_from_seur_lazy)))
    expect_identical(nrow(collect(proximity_from_seur_lazy)), nrow(proximity))
    expect_true(all(colnames(proximity) %in% colnames(proximity_from_seur_lazy)))

    DBI::dbDisconnect(proximity_from_seur_lazy$src$con)
  })

  test_that("ProximityScores fails with invalid input", {
    expect_error(ProximityScores("Invalid"))
    expect_error(ProximityScores(seur_obj, assay = FALSE))
    expect_error(ProximityScores(seur_obj, meta_data_columns = "Invalid"))
    expect_error(ProximityScores(seur_obj, add_marker_counts = "Invalid"))
    expect_error(ProximityScores(seur_obj) <- "Invalid")
  })
}
