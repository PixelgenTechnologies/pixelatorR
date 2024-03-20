options(pixelatorR.arrow_outdir = tempdir())
se <- ReadMPX_Seurat(system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR"),
                     overwrite = TRUE, return_cellgraphassay = TRUE)
options(pixelatorR.interactive = FALSE)

test_that("clean_edgelists_directories works as expected", {
  expect_no_error(clean_edgelists_directories())
})

rm(se)

test_that("clean_edgelists_directories works as expected", {
  expect_no_error(clean_edgelists_directories())
})
