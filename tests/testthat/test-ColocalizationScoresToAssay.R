for (assay_version in c("v3", "v5")) {

  options(Seurat.object.assay.version = assay_version)

  expected_assay_class <- ifelse(assay_version == "v3", "Assay", "Assay5")

  pxl_file <- system.file("extdata/five_cells",
                          "five_cells.pxl",
                          package = "pixelatorR"
  )
  seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

  test_that("ColocalizationScoresToAssay works as expected", {

    # Default settings
    expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj))
    expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization"))
    expect_s4_class(seur_obj[["colocalization"]], expected_assay_class)
    expect_equal(ncol(seur_obj[["colocalization"]]), ncol(seur_obj[["mpxCells"]]))
    expect_equal(
      seur_obj[["colocalization"]]$data[1:2, 1:2],
      new(
        "dgCMatrix",
        i = 0:1,
        p = c(0L, 2L, 2L),
        Dim = c(2L, 2L),
        Dimnames = list(
          c("ACTB-B2M", "ACTB-CD102"),
          c("RCVCMP0000217",
            "RCVCMP0000118")
        ),
        x = c(-0.940132075815013, 1.92413178598625),
        factors = list()
      )
    )

    # New assay name
    expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, new_assay = "coloc"))
    expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization", "coloc"))
    expect_s4_class(seur_obj[["coloc"]], expected_assay_class)
    expect_equal(ncol(seur_obj[["coloc"]]), ncol(seur_obj[["mpxCells"]]))
    expect_equal(nrow(seur_obj[["coloc"]]), 3160)
  })


  test_that("ColocalizationScoresToAssay fails with invalid input", {
    expect_error(seur_obj <- ColocalizationScoresToAssay("Invalid"), "no applicable method")
    expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, assay = "invalid"))
    expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, new_assay = 0), "'new_assay' must be a character of length 1")
    expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, values_from = "invalid"), "'arg' should be one of")
  })

}
