pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expect_no_error(seur_obj <- suppressWarnings(ReadPNA_Seurat(pxl_file, verbose = FALSE)))
  expect_no_error(pna_assay <- seur_obj[["PNA"]])
  expect_no_error(proximity <- ProximityScores(seur_obj))
  expect_no_error(proximity_lazy <- ProximityScores(seur_obj, lazy = TRUE))

  test_that("ProximityScoresToAssay works as expected", {
    # tbl_df
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity))
    expect_s4_class(prox_matrix, "dgCMatrix")
    expect_equal(dim(prox_matrix), c(12561, 5))
    expect_equal(sum(is.na(prox_matrix)), 4109)
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity, missing_obs = 0))
    expect_equal(sum(is.na(prox_matrix)), 0)
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity, return_sparse = FALSE))
    expect_type(prox_matrix, "double")
    expect_true(inherits(prox_matrix, "matrix"))
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity, values_from = "log2_ratio", return_sparse = FALSE))
    expect_equal(mean(prox_matrix, na.rm = T), -0.05025435, tolerance = 1e-6)

    # tbl_lazy
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity, missing_obs = 0))
    expect_no_error(prox_matrix_from_lazy <- ProximityScoresToAssay(proximity_lazy))
    expect_s4_class(prox_matrix_from_lazy, "dgCMatrix")
    expect_equal(dim(prox_matrix_from_lazy), c(12458, 5))
    expect_equal(
      sum(prox_matrix[rownames(prox_matrix_from_lazy), colnames(prox_matrix_from_lazy)] - prox_matrix_from_lazy),
      0
    )

    # PNAAssay
    expect_no_error(prox_assay <- ProximityScoresToAssay(pna_assay))
    expect_s4_class(prox_assay, "Assay5")
    expect_equal(dim(prox_assay), c(12561, 5))

    # Seurat
    expect_no_error(seur_conv <- ProximityScoresToAssay(seur_obj))
    expect_equal(SeuratObject::Assays(seur_conv), c("PNA", "proximity"))
    expect_s4_class(seur_conv, "Seurat")
    DefaultAssay(seur_conv) <- "proximity"
    expect_equal(dim(seur_conv), c(12561, 5))
    expect_equal(
      head(rownames(seur_conv)),
      c(
        "CD56/CD56",
        "CD56/mIgG2b",
        "CD56/CD71",
        "CD56/CD6",
        "CD56/Siglec-9",
        "CD56/CD79a"
      )
    )
  })

  test_that("ProximityScoresToAssay fails with invalid input", {
    # tbl_df
    expect_error(ProximityScoresToAssay("Invalid"))
    expect_error(ProximityScoresToAssay(proximity, values_from = "Invalid"))
    expect_error(ProximityScoresToAssay(proximity, missing_obs = "Invalid"))
    expect_error(ProximityScoresToAssay(proximity, return_sparse = "Invalid"))

    # PNAAssay
    expect_error(ProximityScoresToAssay(prox_assay, values_from = "Invalid"))
    expect_error(ProximityScoresToAssay(prox_assay, missing_obs = "Invalid"))
    expect_error(ProximityScoresToAssay(prox_assay, return_sparse = "Invalid"))

    # Seurat
    expect_error(ProximityScoresToAssay("Invalid"))
    expect_error(ProximityScoresToAssay(seur_obj, values_from = "Invalid"))
    expect_error(ProximityScoresToAssay(seur_obj, missing_obs = "Invalid"))
  })
}
