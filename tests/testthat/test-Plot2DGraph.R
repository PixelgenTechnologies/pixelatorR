pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1:2])
seur_obj[["mpxCells"]] <- KeepLargestComponent(seur_obj[["mpxCells"]])
seur_obj <- ComputeLayout(seur_obj, layout_method = "pmds")

test_that("Plot2DGraph works as expected", {
  expect_no_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds")})
  expect_s3_class(layout_plot, "ggplot")
  expect_no_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", marker = "CD14")})
  expect_equal(layout_plot[[1]]$labels$title, "RCVCMP0000000")
  expect_equal(layout_plot[[1]]$labels$colour, "CD14\n(log-scaled)")
  expect_equal(dim(layout_plot[[1]]$data), c(2613, 8))

  # Test with showBnodes active
  expect_no_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", show_Bnodes = TRUE)})
  expect_equal(dim(layout_plot[[1]]$data), c(3525, 7))

  # Test with return_plot_list = TRUE
  expect_no_error({layout_plots <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1:2], layout_method = "pmds", return_plot_list = TRUE)})
  expect_type(layout_plots, "list")
  expect_equal(length(layout_plots), 2)
})

test_that("Plot2DGraph fails with invalid input", {
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "invalid")})
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", colors = "invalid")},
               "'colors' must be a character vector with at least 2 colors")
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", assay = "invalid")})
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", map_nodes = "invalid")},
               "'map_nodes' must be either TRUE or FALSE")
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", map_edges = "invalid")},
               "'map_edges' must be either TRUE or FALSE")
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", node_size = "invalid")},
               "'node_size' must be a numeric value")
  expect_error({layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", edge_width = "invalid")},
               "'edge_width' must be a numeric value")
})

test_that("Plot2DGraphM works as expected", {
  expect_no_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"))})
  expect_s3_class(layout_plot, "patchwork")
  expect_equal(layout_plot[[1]]$labels$title, NULL)
})

test_that("Plot2DGraph fails with invalid input", {
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "invalid")})
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), colors = "invalid")},
               "'colors' must be a character vector with at least 2 colors")
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), assay = "invalid")})
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), map_nodes = "invalid")},
               "'map_nodes' must be either TRUE or FALSE")
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), map_edges = "invalid")},
               "'map_edges' must be either TRUE or FALSE")
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), node_size = "invalid")},
               "'node_size' must be a numeric value")
  expect_error({layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), edge_width = "invalid")},
               "'edge_width' must be a numeric value")
})