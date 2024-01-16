pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
seur_obj[["mpxCells"]] <- KeepLargestComponent(seur_obj[["mpxCells"]])
seur_obj <- ComputeLayout(seur_obj, layout_method = "pmds")

test_that("Plot2DGraph works as expected", {
  expect_no_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds")})
  expect_s3_class(layout_plot, "ggplot")
  expect_no_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", marker = "CD14")})
})
