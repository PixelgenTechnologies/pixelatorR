options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("PolarizationScoresToAssay works as expected", {

  # Default settings
  expect_no_error(seur_obj <- PolarizationScoresToAssay(seur_obj))
  expect_equal(names(seur_obj@assays), c("mpxCells", "polarization"))
  expect_s4_class(seur_obj[["polarization"]], "Assay")
  expect_equal(dim(seur_obj[["polarization"]]), dim(seur_obj[["mpxCells"]]))
  expect_equal(seur_obj[["polarization"]]$data[1:2, 1:2],
               structure(
                 c(-0.133503001331984,-0.675947218565106, 0, 4.55113839877331),
                 dim = c(2L, 2L),
                 dimnames = list(c("ACTB", "B2M"), c("RCVCMP0000217",
                                                     "RCVCMP0000118"))
               ))

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
  expect_equal(seur_obj[["polarization"]]$data[1:2, 1:2],
               structure(
                 c(
                   -0.00164610432425177,
                   -0.00734402218330873,
                   0,
                   0.0307451831152651
                 ),
                 dim = c(2L, 2L),
                 dimnames = list(c("ACTB", "B2M"), c("RCVCMP0000217",
                                                     "RCVCMP0000118"))
               ))
})


test_that("PolarizationScoresToAssay fails with invalid input", {
  expect_error(seur_obj <- PolarizationScoresToAssay("Invalid"), "no applicable method")
  expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, assay = "invalid"))
  expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, new_assay = 0), "'new_assay' must be a character of length 1")
  expect_error(seur_obj <- PolarizationScoresToAssay(seur_obj, values_from = "invalid"), "'arg' should be one of")
})
