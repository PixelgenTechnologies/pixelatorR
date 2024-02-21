options(pixelatorR.arrow_outdir = tempdir())
pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("ColocalizationScoresToAssay works as expected", {

  # Default settings
  expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj))
  expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization"))
  expect_s4_class(seur_obj[["colocalization"]], "Assay")
  expect_equal(ncol(seur_obj[["colocalization"]]), ncol(seur_obj[["mpxCells"]]))
  exp_matrix_top_left <- matrix(c(6.87611904, -1.02836270, -6.22495456, 0.000000),
                                ncol = 2, byrow = TRUE,
                                dimnames = list(c("ACTB-CD11c", "ACTB-CD14"),
                                                c("RCVCMP0000000", "RCVCMP0000002")))
  expect_equal(exp_matrix_top_left, seur_obj[["colocalization"]]$data[1:2, 1:2])

  # New assay name
  expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, new_assay = "coloc"))
  expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization", "coloc"))
  expect_s4_class(seur_obj[["coloc"]], "Assay")
  expect_equal(ncol(seur_obj[["coloc"]]), ncol(seur_obj[["mpxCells"]]))
  expect_equal(nrow(seur_obj[["coloc"]]), 325)

  # Use jaccard_z values
  expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, values_from = "jaccard_z"))
  expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization", "coloc"))
  expect_s4_class(seur_obj[["colocalization"]], "Assay")
  expect_equal(ncol(seur_obj[["colocalization"]]), ncol(seur_obj[["mpxCells"]]))
  exp_matrix_top_left <- matrix(c(6.9636201306, -0.9352617012, -5.7658771986, 0.0000000000),
                                ncol = 2, byrow = TRUE,
                                dimnames = list(c("ACTB-CD11c", "ACTB-CD14"),
                                                c("RCVCMP0000000", "RCVCMP0000002")))
  expect_equal(exp_matrix_top_left, seur_obj[["colocalization"]]$data[1:2, 1:2])
})


test_that("ColocalizationScoresToAssay fails with invalid input", {
  expect_error(seur_obj <- ColocalizationScoresToAssay("Invalid"), "no applicable method")
  expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, assay = "invalid"))
  expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, new_assay = 0), "'new_assay' must be a character of length 1")
  expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, values_from = "invalid"), "'arg' should be one of")
})
