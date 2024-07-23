for (assay_version in c("v3", "v5")) {

  options(Seurat.object.assay.version = assay_version)

  pxl_file <- system.file("extdata/five_cells",
                          "five_cells.pxl",
                          package = "pixelatorR")
  seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
  seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1:2])
  seur_obj <- ComputeLayout(seur_obj, layout_method = "pmds", dim = 3)

  test_that("Plot3DGraph works as expected", {
    expect_no_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", marker = "CD14")})
    expect_no_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", marker = "CD14")})
    expect_s3_class(layout_plot, "plotly")
    expect_no_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", marker = "CD14")})
    expect_equal(layout_plot$x$layoutAttrs[[1]]$annotations$text, "CD14")

    # Test with showBnodes active
    expect_no_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", show_Bnodes = TRUE, marker = "CD14")})

    # Test with project active
    expect_no_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", project = TRUE, marker = "CD14")})
  })

  test_that("Plot3DGraph fails with invalid input", {
    expect_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "invalid", marker = "CD14")})
    expect_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", colors = c("red"), marker = "CD14")},
                 "'colors' must be a character vector with at least 2 color names")
    expect_error({layout_plot <- Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1:2], layout_method = "pmds_3d", node_size = 2, marker = "CD14")},
                 "'cell_id' must be a non-empty character vector with a single cell ID")
  })

}
