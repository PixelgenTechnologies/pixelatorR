# Load example data as a Seurat object
pxl_file_mpx <- minimal_mpx_pxl_file()
seur_obj_mpx <- ReadMPX_Seurat(pxl_file_mpx)
pxl_file_pna <- minimal_pna_pxl_file()
seur_obj_pna <- ReadPNA_Seurat(pxl_file_pna)

test_that("MoleculeRankPlot works for Seurat objects", {
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_mpx)
  })
  expect_s3_class(moleculerank_plot, "ggplot")
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_mpx, group_by = "leiden")
  })
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_pna)
  })
  expect_s3_class(moleculerank_plot, "ggplot")
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_pna, group_by = "tau_type")
  })

  # min max threshold
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_pna, n_umi_min_threshold = 100, n_umi_max_threshold = 1000)
  })

  # facet
  seur_obj_pna$sample_id <- c("S1", "S2", "S2", "S2", "S2")
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_pna, group_by = "sample_id", split = TRUE)
  })
  expect_equal(moleculerank_plot@facet$params$facets %>% names(), "sample_id")

  # rug
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_pna, rug = TRUE)
  })

  # highlight cell counts
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_pna, highlight_cell_counts = TRUE)
  })
})

test_that("MoleculeRankPlot works for data.frame-like objects", {
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj_mpx[[]])
  })
  expect_s3_class(moleculerank_plot, "ggplot")
})
