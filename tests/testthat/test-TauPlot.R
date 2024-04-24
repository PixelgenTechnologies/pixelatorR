pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file)

test_that("TauPlot works for Seurat objects", {
  expect_no_error({tau_plot <- TauPlot(seur_obj)})
  expect_s3_class(tau_plot, "ggplot")
})

test_that("TauPlot works for data.frame-like objects", {
  expect_no_error({tau_plot <- TauPlot(seur_obj[[]])})
  expect_s3_class(tau_plot, "ggplot")
})

test_that("TauPlot fails with invalid input", {
  expect_error(TauPlot("invalid"))
  expect_error(TauPlot(seur_obj, group_by = "invalid"))
})
