options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/mock_data",
                        "mock_mpx_data.pxl",
                        package = "pixelatorR")

seur_obj <- ReadMPX_Seurat(pxl_file, return_cellgraphassay = TRUE, overwrite = TRUE)
outfile <- tempfile()

test_that("WriteMPX.CellGraphAssay works as expected", {
  expect_no_error(WriteMPX(object = seur_obj[["mpxCells"]], outfile, overwrite = TRUE))
  se_obj <- readRDS(outfile)
  expect_no_error(WriteMPX(object = seur_obj[["mpxCells"]], outfile, overwrite = TRUE))
})

test_that("WriteMPX.CellGraphAssay fails when expected", {
  expect_error({suppressWarnings({WriteMPX(object = se_obj, "Invalid/path", overwrite = TRUE)})})
  expect_error({suppressWarnings({WriteMPX(object = se_obj, NULL)})})
})

test_that("WriteMPX.Seurat works as expected", {
  expect_no_error(WriteMPX(object = seur_obj, outfile, overwrite = TRUE))
})

test_that("WriteMPX.Seurat fails when expected", {
  expect_error({suppressWarnings({WriteMPX(object = seur_obj, "Invalid/path", overwrite = TRUE)})})
  expect_error({suppressWarnings({WriteMPX(object = seur_obj, NULL)})})
})
