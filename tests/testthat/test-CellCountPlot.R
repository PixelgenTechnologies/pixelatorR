options(pixelatorR.arrow_outdir = tempdir())
pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("CellCountPlot works for Seurat objects", {
  expect_no_error({cellcount_plot <- CellCountPlot(seur_obj, color_by = "leiden")})
  expect_s3_class(cellcount_plot, "ggplot")
  expect_no_error({cellcount_plot <- CellCountPlot(seur_obj, group_by = "leiden", color_by = "tau_type")})
})

test_that("CellCountPlot works for data.frame-like objects", {
  expect_no_error({cellcount_plot <- CellCountPlot(seur_obj[[]], color_by = "leiden")})
  expect_s3_class(cellcount_plot, "ggplot")
})

test_that("CellCountPlot fails with invalid input", {
  expect_error(CellCountPlot(seur_obj), 'argument "color_by" is missing, with no default')
  expect_error(CellCountPlot(seur_obj, color_by = "edges"), "'color_by' must be a character or factor")
  expect_error(CellCountPlot(seur_obj, color_by = "leiden", group_by = "edges"), "'group_by' must be a character or factor")
})
