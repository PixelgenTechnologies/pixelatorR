pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expect_no_error(seur_obj <- suppressWarnings(ReadPNA_Seurat(pxl_file, verbose = FALSE)))
  expect_no_error(pna_assay <- seur_obj[["PNA"]])
  expect_no_error(proximity <-
    ProximityScores(seur_obj) %>%
    filter(log2_ratio != 0))
  expect_no_error(proximity_lazy <-
    ProximityScores(seur_obj, lazy = TRUE) %>%
    filter(log2_ratio != 0))

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
        "CD56:HLA-DR",
        "CD56:CD9",
        "CD56:CD59"
      )
    )

    # Test separator
    expect_no_error(prox_matrix <- ProximityScoresToAssay(proximity, separator = "_"))

    # Test that all-zero components are included with value 0
    proximity_scores <-
      structure(
        list(
          marker_1 = c(
            "CD2", "CD2", "CD2", "CD2", "CD2",
            "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2",
            "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2",
            "CD2", "CD2", "CD2"
          ),
          marker_2 = c(
            "CD2", "CD3e", "CD2", "CD3e",
            "CD4", "CD2", "CD3e", "CD2", "CD2", "CD3e", "CD2", "CD4", "CD3e",
            "CD4", "CD2", "CD2", "CD3e", "CD2", "CD2", "CD4", "CD3e", "CD2",
            "CD3e", "CD4", "CD2", "CD2"
          ),
          component = c(
            "S1_7e9a3f5ef0362b6d",
            "S1_7e9a3f5ef0362b6d", "S1_99aa508551f33cd3", "S1_99aa508551f33cd3",
            "S1_99aa508551f33cd3", "S1_2ee5115016f6d3bf", "S1_2ee5115016f6d3bf",
            "S1_dc8d14e01732603e", "S1_c8d43c718f2181fc", "S1_c8d43c718f2181fc",
            "S1_844fa3a723a26b8a", "S1_844fa3a723a26b8a", "S1_a0525046264c9183",
            "S1_a0525046264c9183", "S1_a0525046264c9183", "S1_128b348aca07cb57",
            "S1_128b348aca07cb57", "S1_e2055911bf1693bd", "S1_f0dceb55e0820e2d",
            "S1_f0dceb55e0820e2d", "S1_f0dceb55e0820e2d", "S1_b720a4c4f4bfac3e",
            "S1_b720a4c4f4bfac3e", "S1_b720a4c4f4bfac3e", "S1_d375cd40d5330b42",
            "S1_a9fec42b5a90b80a"
          ),
          log2_ratio = c(
            0, 0, 0.84546704897467,
            -0.0768522364596248, 0.223438231487999, 0, 0, 0, 0.24297675349254,
            0.0657980509150873, 0, 0, -1.05763009717606, 2.09541956507868,
            0, 0.761213140412883, -0.203913000728787, 0, -0.50927356912461,
            0.362157939675895, -1.53439367631852, 1, -0.193546480683928,
            0.305143510591799, 0, 0
          )
        ),
        row.names = c(NA, -26L),
        class = c(
          "tbl_df",
          "tbl", "data.frame"
        )
      )

    expect_warning(
      proximity_scores_wide <- ProximityScoresToAssay(
        proximity_scores,
        values_from = "log2_ratio"
      ),
      "only zero proximity scores"
    )

    expect_true(
      all(unique(proximity_scores$component) %in% colnames(proximity_scores_wide))
    )
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
