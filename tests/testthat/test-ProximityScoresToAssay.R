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
    expected_result <- new(
      "dgCMatrix",
      i = c(0L, 1L, 1L),
      p = c(0L, 2L, 2L, 2L, 2L, 3L),
      Dim = c(2L, 5L),
      Dimnames = list(
        c("CD56:CD6", "CD56:HLA-ABC"),
        c(
          "c3c393e9a17c1981",
          "d4074c845bb62800",
          "efe0ed189cb499fc",
          "0a45497c6bfbfb22",
          "2708240b908e2eba"
        )
      ),
      x = c(-1.0496307677246, -1.66902676550963, -1.13093086982645),
      factors = list()
    )
    expect_equal(prox_matrix %>% head(2), expected_result)

    # tbl_lazy
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity))
    expect_no_error(prox_matrix_from_lazy <- ProximityScoresToAssay(proximity_lazy))
    expect_s4_class(prox_matrix_from_lazy, "dgCMatrix")
    expect_equal(dim(prox_matrix_from_lazy), c(6057, 5))
    expect_equal(
      dim(prox_matrix), dim(prox_matrix_from_lazy)
    )
    expect_equal(sum(prox_matrix) - sum(prox_matrix_from_lazy), 0) # should be 0

    # PNAAssay
    expect_no_error(prox_assay <- ProximityScoresToAssay(pna_assay))
    expect_s4_class(prox_assay, "Assay5")
    expect_equal(dim(prox_assay), c(6057, 5))

    # Seurat
    expect_no_error(seur_conv <- ProximityScoresToAssay(seur_obj))
    expect_equal(SeuratObject::Assays(seur_conv), c("PNA", "proximity"))
    expect_s4_class(seur_conv, "Seurat")
    DefaultAssay(seur_conv) <- "proximity"
    expect_equal(dim(seur_conv), c(6057, 5))
    expect_equal(
      head(rownames(seur_conv)),
      c(
        "CD56:CD6",
        "CD56:HLA-ABC",
        "CD56:CD80",
        "CD56:CD9",
        "CD56:CD59",
        "CD56:HLA-DR-DP-DQ"
      )
    )

    # Test separator
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity, separator = "_"))
  })

  test_that("ProximityScoresToAssay fails with invalid input", {
    # tbl_df
    expect_error(ProximityScoresToAssay("Invalid"))
    expect_error(ProximityScoresToAssay(proximity, values_from = "Invalid"))

    # PNAAssay
    expect_error(ProximityScoresToAssay(prox_assay, values_from = "Invalid"))

    # Seurat
    expect_error(ProximityScoresToAssay("Invalid"))
    expect_error(ProximityScoresToAssay(seur_obj, values_from = "Invalid"))

    # separator must be single character
    expect_error(ProximityScoresToAssay(proximity, separator = "Invalid"))

    # separator must be absent from marker names
    expect_error(ProximityScoresToAssay(proximity, separator = "-"))
  })
}
