library(pixelatorR)
pxl_file <- system.file("extdata/PBMC_10_cells",
  "Sample01_test.pxl",
  package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("PolarizationScoresToAssay works as expected", {

  # Default settings
  expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj))
  expect_equal(names(seur_obj@assays), c("mpxCells", "polarization"))
  expect_s4_class(seur_obj[["polarization"]], "Assay")
  expect_equal(dim(seur_obj[["polarization"]]), dim(seur_obj[["mpxCells"]]))
  exp_matrix_top_left <- matrix(c(-0.02884569, 0.01942594, 0.64174516, 2.81220649),
                       ncol = 2, byrow = TRUE,
                       dimnames = list(c("ACTB", "CD11c"),
                                       c("RCVCMP0000000", "RCVCMP0000002")))
  expect_equal(exp_matrix_top_left, seur_obj[["polarization"]]$data[1:2, 1:2])

  # New assay name
  expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj, new_assay = "pol"))
  expect_equal(names(seur_obj@assays), c("mpxCells", "polarization", "pol"))
  expect_s4_class(seur_obj[["pol"]], "Assay")
  expect_equal(dim(seur_obj[["pol"]]), dim(seur_obj[["mpxCells"]]))

  # Use morans_i values
  expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj, values_from = "morans_i"))
  expect_equal(names(seur_obj@assays), c("mpxCells", "polarization", "pol"))
  expect_s4_class(seur_obj[["polarization"]], "Assay")
  expect_equal(dim(seur_obj[["polarization"]]), dim(seur_obj[["mpxCells"]]))
  exp_matrix_top_left <- matrix(c(-0.0004683552, -0.0005655225, 0.0023255354, 0.0172761591),
                                ncol = 2, byrow = TRUE,
                                dimnames = list(c("ACTB", "CD11c"),
                                                c("RCVCMP0000000", "RCVCMP0000002")))
  expect_equal(exp_matrix_top_left, seur_obj[["polarization"]]$data[1:2, 1:2])
})


test_that("PolarizationScoresToAssay fails with invalid input", {
  expect_error(seur_obj <- PolarizationScoresToAssay("Invalid"), "no applicable method")
  expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, assay = "invalid"))
  expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, new_assay = 0), "'new_assay' must be a character of length 1")
  expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, values_from = "invalid"), "'arg' should be one of")
})
