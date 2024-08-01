for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expected_assay_class <- ifelse(assay_version == "v3", "Assay", "Assay5")

  pxl_file <- system.file("extdata/five_cells",
    "five_cells.pxl",
    package = "pixelatorR"
  )
  seur_obj <- ReadMPX_Seurat(pxl_file)

  test_that("PolarizationScoresToAssay works as expected", {
    # Default settings
    expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj))
    expect_equal(names(seur_obj@assays), c("mpxCells", "polarization"))
    expect_s4_class(seur_obj[["polarization"]], expected_assay_class)
    expect_equal(dim(seur_obj[["polarization"]]), dim(seur_obj[["mpxCells"]]))
    expect_equal(
      seur_obj[["polarization"]]$data[1:2, 1:2],
      new(
        "dgCMatrix",
        i = c(0L, 1L, 1L),
        p = c(0L, 2L, 3L),
        Dim = c(
          2L,
          2L
        ),
        Dimnames = list(c("ACTB", "B2M"), c("RCVCMP0000217", "RCVCMP0000118")),
        x = c(-0.133503001331984, -0.675947218565106, 4.55113839877331),
        factors = list()
      )
    )

    # New assay name
    expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj, new_assay = "pol"))
    expect_equal(names(seur_obj@assays), c("mpxCells", "polarization", "pol"))
    expect_s4_class(seur_obj[["pol"]], expected_assay_class)
    expect_equal(dim(seur_obj[["pol"]]), dim(seur_obj[["mpxCells"]]))

    # Use morans_i values
    expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj, values_from = "morans_i"))
    expect_equal(names(seur_obj@assays), c("mpxCells", "polarization", "pol"))
    expect_s4_class(seur_obj[["polarization"]], expected_assay_class)
    expect_equal(dim(seur_obj[["polarization"]]), dim(seur_obj[["mpxCells"]]))
    expect_equal(
      seur_obj[["polarization"]]$data[1:2, 1:2],
      new(
        "dgCMatrix",
        i = c(0L, 1L, 1L),
        p = c(0L, 2L, 3L),
        Dim = c(
          2L,
          2L
        ),
        Dimnames = list(c("ACTB", "B2M"), c("RCVCMP0000217", "RCVCMP0000118")),
        x = c(
          -0.00164610432425177,
          -0.00734402218330873,
          0.0307451831152651
        ),
        factors = list()
      )
    )

    # Use dashes in component IDs
    expect_no_error({
      seur_obj <- SeuratObject::RenameCells(seur_obj, new.names = paste0("A-1_", colnames(seur_obj)))
    })
    expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj))
  })


  test_that("PolarizationScoresToAssay fails with invalid input", {
    expect_error(seur_obj <- PolarizationScoresToAssay("Invalid"), "no applicable method")
    expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, assay = "invalid"))
    expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, new_assay = 0), "'new_assay' must be a character of length 1")
    expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, values_from = "invalid"), "'arg' should be one of")
  })
}
