pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("CellCountPlot works for Seurat objects", {
  expect_no_error({cellcount_plot <- CellCountPlot(seur_obj, color.by = "leiden")})
  expect_s3_class(cellcount_plot, "ggplot")
  expect_no_error({cellcount_plot <- CellCountPlot(seur_obj, group.by = "leiden", color.by = "tau_type")})
})

test_that("CellCountPlot works for data.frame-like objects", {
  expect_no_error({cellcount_plot <- CellCountPlot(seur_obj[[]], color.by = "leiden")})
  expect_s3_class(cellcount_plot, "ggplot")
})

test_that("CellCountPlot fails with invalid input", {
  expect_error(CellCountPlot(seur_obj), 'argument "color.by" is missing, with no default')
  expect_error(CellCountPlot(seur_obj, color.by = "edges"), "'color.by' must be a character or factor")
  expect_error(CellCountPlot(seur_obj, color.by = "leiden", group.by = "edges"), "'group.by' must be a character or factor")
})
