options(pixelatorR.arrow_outdir = tempdir())
pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("ColocalizationScoresToAssay works as expected", {

  # Default settings
  expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj))
  expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization"))
  expect_s4_class(seur_obj[["colocalization"]], "Assay")
  expect_equal(ncol(seur_obj[["colocalization"]]), ncol(seur_obj[["mpxCells"]]))
  expect_equal(seur_obj[["colocalization"]]$data[1:2, 1:2],
               structure(
                 c(-0.940132075815013, 1.92413178598625, 0, 0),
                 dim = c(2L,
                         2L),
                 dimnames = list(
                   c("ACTB-B2M", "ACTB-CD102"),
                   c("RCVCMP0000217",
                     "RCVCMP0000118")
                 )
               ))

  # New assay name
  expect_no_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, new_assay = "coloc"))
  expect_equal(names(seur_obj@assays), c("mpxCells", "colocalization", "coloc"))
  expect_s4_class(seur_obj[["coloc"]], "Assay")
  expect_equal(ncol(seur_obj[["coloc"]]), ncol(seur_obj[["mpxCells"]]))
  expect_equal(nrow(seur_obj[["coloc"]]), 3160)
})


test_that("ColocalizationScoresToAssay fails with invalid input", {
  expect_error(seur_obj <- ColocalizationScoresToAssay("Invalid"), "no applicable method")
  expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, assay = "invalid"))
  expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, new_assay = 0), "'new_assay' must be a character of length 1")
  expect_error(seur_obj <- ColocalizationScoresToAssay(seur_obj, values_from = "invalid"), "'arg' should be one of")
})
