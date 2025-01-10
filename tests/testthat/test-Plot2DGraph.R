for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  pxl_file <- system.file("extdata/five_cells",
    "five_cells.pxl",
    package = "pixelatorR"
  )
  seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
  seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1:2])
  seur_obj <- ComputeLayout(seur_obj, layout_method = "pmds")

  test_that("Plot2DGraph works as expected", {
    expect_no_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds")
    })
    expect_s3_class(layout_plot, "ggplot")
    expect_no_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", marker = "CD3E")
    })
    expect_equal(layout_plot[[1]]$labels$title, "RCVCMP0000217")
    expect_equal(layout_plot[[1]]$labels$colour, "CD3E\n(log-scaled)")
    expect_equal(dim(layout_plot[[1]]$data), c(1395, 8))

    # Test with showBnodes active
    expect_no_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", show_Bnodes = TRUE)
    })
    expect_equal(dim(layout_plot[[1]]$data), c(2470, 7))

    # Test with return_plot_list = TRUE
    expect_no_error({
      layout_plots <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1:2], layout_method = "pmds", return_plot_list = TRUE)
    })
    expect_type(layout_plots, "list")
    expect_equal(length(layout_plots), 2)
  })

  test_that("Plot2DGraph fails with invalid input", {
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", colors = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", assay = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", map_nodes = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", map_edges = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", node_size = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraph(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", edge_width = "invalid")
    })
  })

  test_that("Plot2DGraphM works as expected", {
    expect_no_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"))
    })
    expect_s3_class(layout_plot, "patchwork")
    expect_equal(layout_plot[[1]]$labels$title, NULL)
  })

  test_that("Plot2DGraph fails with invalid input", {
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), colors = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), assay = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), map_nodes = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), map_edges = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), node_size = "invalid")
    })
    expect_error({
      layout_plot <- Plot2DGraphM(seur_obj, cells = colnames(seur_obj)[1], layout_method = "pmds", markers = c("HLA-DR"), edge_width = "invalid")
    })
  })
}
