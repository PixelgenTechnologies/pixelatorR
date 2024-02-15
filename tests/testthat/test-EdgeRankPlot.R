# Load example data as a Seurat object
pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- suppressWarnings({ReadMPX_Seurat(pxl_file, overwrite = TRUE)})

test_that("EdgeRankPlot works for Seurat objects", {
  expect_no_error({edgerank_plot <- EdgeRankPlot(seur_obj)})
  expect_s3_class(edgerank_plot, "ggplot")
  expect_no_error({edgerank_plot <- EdgeRankPlot(seur_obj, group_by = "leiden")})
})

test_that("EdgeRankPlot works for data.frame-like objects", {
  expect_no_error({edgerank_plot <- EdgeRankPlot(seur_obj[[]])})
  expect_s3_class(edgerank_plot, "ggplot")
})

