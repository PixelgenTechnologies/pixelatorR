options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/mock_data",
                        "mock_mpx_data.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

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
