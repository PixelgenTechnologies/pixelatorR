library(Seurat)
pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)

seur_obj <- NormalizeData(seur_obj, normalization.method = "LogNormalize")



test_that("AbundanceColocalizationPlot works as expected", {

  expect_no_error(AbundanceColocalizationPlot(seur_obj,
                                              markers_x = "CD20",
                                              markers_y = "CD3E",
                                              layer = "counts"))


  expect_no_error(AbundanceColocalizationPlot(seur_obj,
                              markers_x = c("CD20", "CD19", "CD22"),
                              markers_y = c("CD2", "CD3E", "CD5"),
                              shared_scales = FALSE,
                              coord_fixed = FALSE,
                              draw_origo = FALSE,
                              coloc_score = "pearson",
                              layer = "counts",
                              pt_size = c(0.1, 2)))

  expect_error(AbundanceColocalizationPlot(seur_obj,
                                              markers_x = c("CD20", "MarkerNotHere", "CD22"),
                                              markers_y = c("CD2", "CD3E", "CD5")))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19", "CD22"),
                                           markers_y = c("CD2", "MarkerNotHere", "CD5")))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19"),
                                           markers_y = c("CD2", "CD5"),
                                           shared_scales = "FALSE"))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19"),
                                           markers_y = c("CD2", "CD5"),
                                           coord_fixed = "FALSE"))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19"),
                                           markers_y = c("CD2", "CD5"),
                                           draw_origo = "FALSE"))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19"),
                                           markers_y = c("CD2", "CD5"),
                                           coloc_score = "ColocScoreNotHere"))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19"),
                                           markers_y = c("CD2", "CD5"),
                                           layer = "LayerNotHere"))
  expect_error(AbundanceColocalizationPlot(seur_obj,
                                           markers_x = c("CD20", "CD19"),
                                           markers_y = c("CD2", "CD5"),
                                           pt_size = c(1, 2, 3)))
})
