options(pixelatorR.arrow_outdir = tempdir())
se <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
                     overwrite = TRUE, return_cellgraphassay = TRUE)
options(pixelatorR.interactive = FALSE)

test_that("edgelist_directories_clean works as expected", {
  expect_no_error(edgelist_directories_clean())
})

rm(se)

test_that("edgelist_directories_clean works as expected", {
  expect_no_error(edgelist_directories_clean())
})
