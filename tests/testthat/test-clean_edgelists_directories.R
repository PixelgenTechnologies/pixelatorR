options(pixelatorR.arrow_outdir = tempdir())
options(pixelatorR.interactive = FALSE)
rm(list = ls())
edgelist_directories_clean()

se <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
                     overwrite = TRUE, return_cellgraphassay = TRUE)

test_that("edgelist_directories_clean works as expected", {
  expect_no_error(edgelist_directories_clean())
})

rm(se)
test_that("edgelist_directories_clean works as expected", {
  expect_no_error(edgelist_directories_clean())
})

se <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
                     overwrite = TRUE, return_cellgraphassay = TRUE)

test_that("edgelist_directories_du works as expected", {
  expect_no_error(du <- edgelist_directories_du())
  expect_equal(du, fs::fs_bytes(487466))
  expect_no_error(du_all <- edgelist_directories_du(list_all = TRUE))
})
