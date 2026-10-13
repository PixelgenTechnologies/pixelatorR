library(Seurat)
pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  seur_obj <- ReadMPX_Seurat(pxl_file)

  seur_obj <- NormalizeData(seur_obj, normalization.method = "LogNormalize")

  if (assay_version == "v3") {
    test_that("AbundanceColocalizationPlot warns once per session that it will be removed", {
      rlang::local_options(lifecycle_verbosity = "default")

      expect_warning(
        expect_error(AbundanceColocalizationPlot(
          seur_obj,
          markers_x = "MarkerNotHere",
          markers_y = "CD3E"
        )),
        regexp = "will be removed in a future release",
        class = "lifecycle_warning_deprecated"
      )

      expect_no_condition(
        plot <- AbundanceColocalizationPlot(
          seur_obj,
          markers_x = "CD20",
          markers_y = "CD3E",
          layer = "counts"
        ),
        class = "lifecycle_warning_deprecated"
      )
      expect_s3_class(plot, "ggplot")
    })
  }

  test_that("AbundanceColocalizationPlot works as expected", {
    expect_no_error(lifecycle::expect_deprecated(AbundanceColocalizationPlot(seur_obj,
      markers_x = "CD20",
      markers_y = "CD3E",
      layer = "counts"
    )))


    expect_no_error(lifecycle::expect_deprecated(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19", "CD22"),
      markers_y = c("CD2", "CD3E", "CD5"),
      shared_scales = FALSE,
      coord_fixed = FALSE,
      draw_origo = FALSE,
      coloc_score = "pearson",
      layer = "counts",
      pt_size = c(0.1, 2)
    )))

    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "MarkerNotHere", "CD22"),
      markers_y = c("CD2", "CD3E", "CD5")
    )))
    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19", "CD22"),
      markers_y = c("CD2", "MarkerNotHere", "CD5")
    )))
    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19"),
      markers_y = c("CD2", "CD5"),
      shared_scales = "FALSE"
    )))
    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19"),
      markers_y = c("CD2", "CD5"),
      coord_fixed = "FALSE"
    )))
    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19"),
      markers_y = c("CD2", "CD5"),
      draw_origo = "FALSE"
    )))
    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19"),
      markers_y = c("CD2", "CD5"),
      coloc_score = "ColocScoreNotHere"
    )))

    lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
      markers_x = c("CD20", "CD19"),
      markers_y = c("CD2", "CD5"),
      pt_size = c(1, 2, 3)
    )))

    if (assay_version == "v5") {
      # Seurat also warns about the missing layer. Quiet the deprecation
      # signal so this expectation still checks that layer warning.
      rlang::local_options(lifecycle_verbosity = "quiet")
      expect_warning(
        expect_error(AbundanceColocalizationPlot(seur_obj,
          markers_x = c("CD20", "CD19"),
          markers_y = c("CD2", "CD5"),
          layer = "LayerNotHere"
        ))
      )
    } else {
      lifecycle::expect_deprecated(expect_error(AbundanceColocalizationPlot(seur_obj,
        markers_x = c("CD20", "CD19"),
        markers_y = c("CD2", "CD5"),
        layer = "LayerNotHere"
      )))
    }
  })
}
