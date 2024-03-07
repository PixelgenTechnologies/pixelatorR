options(pixelatorR.arrow_outdir = tempdir())

test_that("CellGraphs.Seurat getter/setter works as expected", {
  se <- ReadMPX_Seurat(system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR"),
                                   overwrite = TRUE, return_cellgraphassay = TRUE)
  cg_list <- CellGraphs(se)
  expect_type(cg_list, "list")
  expect_equal(cg_list %>% length(), 5)
  CellGraphs(se) <- cg_list
})

test_that("CellGraphs.Seurat getter/setter fails when invalid input is provided", {
  se <- ReadMPX_Seurat(system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR"),
                                   overwrite = TRUE, return_cellgraphassay = TRUE)
  expect_error(CellGraphs("Invalid input"))
  expect_error(CellGraphs(se) <- "Invalid input", "Invalid class 'character'")
})
