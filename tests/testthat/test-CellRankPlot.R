# Load example data as a Seurat object
pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("CellRankPlot works for Seurat objects", {
  expect_no_error({edgerank_plot <- CellRankPlot(seur_obj)})
  expect_s3_class(edgerank_plot, "ggplot")
  expect_no_error({edgerank_plot <- CellRankPlot(seur_obj, group_by = "leiden")})
})

test_that("CellRankPlot works for data.frame-like objects", {
  expect_no_error({edgerank_plot <- CellRankPlot(seur_obj[[]])})
  expect_s3_class(edgerank_plot, "ggplot")
})

