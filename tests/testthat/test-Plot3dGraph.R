for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  pxl_file <- minimal_mpx_pxl_file()
  seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
  seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1:2])
  seur_obj <- ComputeLayout(seur_obj, layout_method = "pmds")

  test_that("Plot3DGraph works as expected", {
    layout_plot <- lifecycle::expect_deprecated(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", marker = "CD14")
    )
    layout_plot <- lifecycle::expect_deprecated(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", marker = "CD14")
    )
    expect_s3_class(layout_plot, "plotly")
    layout_plot <- lifecycle::expect_deprecated(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", marker = "CD14")
    )
    expect_equal(layout_plot$x$layoutAttrs[[1]]$annotations$text, "CD14")

    # Test with showBnodes active
    lifecycle::expect_deprecated(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", show_Bnodes = TRUE, marker = "CD14")
    )

    # Test with project active
    lifecycle::expect_deprecated(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", project = TRUE, marker = "CD14")
    )
  })

  test_that("Plot3DGraph fails with invalid input", {
    lifecycle::expect_deprecated(expect_error(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "invalid", marker = "CD14")
    ))
    lifecycle::expect_deprecated(expect_error(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1], layout_method = "pmds_3d", colors = c("red"), marker = "CD14")
    ))
    lifecycle::expect_deprecated(expect_error(
      Plot3DGraph(seur_obj, cell_id = colnames(seur_obj)[1:2], layout_method = "pmds_3d", node_size = 2, marker = "CD14")
    ))
  })
}
