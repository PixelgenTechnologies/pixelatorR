pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

test_that("pxContent_vs_Tau works for Seurat objects", {
  expect_no_error({tau_plot <- pxContent_vs_Tau(seur_obj)})
  expect_s3_class(tau_plot, "ggplot")
})

test_that("pxContent_vs_Tau works for data.frame-like objects", {
  expect_no_error({tau_plot <- pxContent_vs_Tau(seur_obj[[]])})
  expect_s3_class(tau_plot, "ggplot")
})

test_that("pxContent_vs_Tau fails with invalid input", {
  expect_error(pxContent_vs_Tau("invalid"))
  expect_error(pxContent_vs_Tau(seur_obj, group.by = "invalid"))
})
