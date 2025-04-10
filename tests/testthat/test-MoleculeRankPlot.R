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
})

test_that("MoleculeRankPlot works for data.frame-like objects", {
  expect_no_error({
    moleculerank_plot <- MoleculeRankPlot(seur_obj[[]])
  })
  expect_s3_class(moleculerank_plot, "ggplot")
})
