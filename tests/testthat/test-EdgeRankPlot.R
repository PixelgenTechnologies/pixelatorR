options(pixelatorR.arrow_outdir = tempdir())

# Load example data as a Seurat object
pxl_file <- system.file("extdata/mock_data",
                        "mock_mpx_data.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("EdgeRankPlot works for Seurat objects", {
  expect_no_error({edgerank_plot <- EdgeRankPlot(seur_obj)})
  expect_s3_class(edgerank_plot, "ggplot")
  expect_no_error({edgerank_plot <- EdgeRankPlot(seur_obj, group_by = "leiden")})
})

test_that("EdgeRankPlot works for data.frame-like objects", {
  expect_no_error({edgerank_plot <- EdgeRankPlot(seur_obj[[]])})
  expect_s3_class(edgerank_plot, "ggplot")
})

