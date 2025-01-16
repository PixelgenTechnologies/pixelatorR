library(dplyr)

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  se <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
    overwrite = TRUE, return_cellgraphassay = TRUE
  )

  se <- LoadCellGraphs(se, cells = colnames(se)[1:2])

  test_that("ComputeLayout works as expected", {
    # Seurat object
    expect_no_error(se <- se %>% ComputeLayout())

    # CellGraphAssay
    layout_method <- "pmds"
    expect_no_error(cg_assay <- se[["mpxCells"]] %>% ComputeLayout(layout_method = layout_method))
    expect_true(layout_method %in% names(cg_assay@cellgraphs[[1]]@layout))
    expect_equal(c("x", "y"), colnames(cg_assay@cellgraphs[[1]]@layout[[layout_method]]))

    # Test with three dimensions
    layout_method_3d <- "pmds_3d"
    expect_no_error(cg_assay <- se[["mpxCells"]] %>%
      ComputeLayout(layout_method = layout_method, dim = 3))
    expect_equal(c("x", "y", "z"), colnames(cg_assay@cellgraphs[[1]]@layout[[layout_method_3d]]))

    # Test with normalize_layout
    expect_no_error(cg_assay <- se[["mpxCells"]] %>%
      ComputeLayout(dim = 3, layout_method = layout_method, normalize_layout = TRUE))
    median_radius <- cg_assay@cellgraphs[[1]]@layout[[layout_method_3d]] %>%
      mutate(across(x:z, ~ .x^2)) %>%
      rowSums() %>%
      sqrt() %>%
      median()
    expect_equal(median_radius, 1)

    # Test with project_on_unit_sphere
    expect_no_error(cg_assay <- se[["mpxCells"]] %>%
      ComputeLayout(dim = 3, layout_method = layout_method, project_on_unit_sphere = TRUE))
    radii <- cg_assay@cellgraphs[[1]]@layout[[layout_method_3d]] %>%
      mutate(across(x:z, ~ .x^2)) %>%
      rowSums() %>%
      sqrt() %>%
      max()
    expect_equal(unique(radii) %>% length(), 1)

    # Weighted pmds
    layout_method <- "wpmds"
    layout_method_3d <- "wpmds_3d"
    expect_no_error(se <- se %>% ComputeLayout(layout_method = layout_method, dim = 3))
    expect_true(layout_method_3d %in% names(CellGraphs(se)[[1]]@layout))
    expect_equal(c("x", "y", "z"), colnames(CellGraphs(se)[[1]]@layout[[layout_method_3d]]))
  })

  test_that("ComputeLayout works as expected with a custom layout function", {
    custom_layout_fkn <- graphlayouts::layout_with_pmds

    cg <- CellGraphs(se)[[colnames(se)[1]]]

    # tbl_graph
    expect_no_error(layout <- ComputeLayout(cg@cellgraph, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
    expect_equal(dim(layout), c(2470, 2))

    # CellGraph
    expect_no_error(cg_layout <- ComputeLayout(cg, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
    expect_true("custom" %in% names(cg_layout@layout))
    expect_equal(dim(cg_layout@layout[["custom"]]), c(2470, 2))

    # CellGraphAssay
    expect_no_error(cg_assay_layout <- ComputeLayout(se[["mpxCells"]], custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
    expect_true("custom" %in% names(CellGraphs(cg_assay_layout)[[1]]@layout))
    expect_equal(dim(CellGraphs(cg_assay_layout)[[1]]@layout[["custom"]]), c(2470, 2))

    # Seurat
    expect_no_error(se_layout <- ComputeLayout(se, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
    expect_true("custom" %in% names(CellGraphs(se_layout)[[1]]@layout))
    expect_equal(dim(CellGraphs(se_layout)[[1]]@layout[["custom"]]), c(2470, 2))

    # Test with new layout name
    expect_no_error(se_layout <- ComputeLayout(se,
      custom_layout_function = custom_layout_fkn,
      custom_layout_function_args = list(pivots = 100), layout_name = "my_layout"
    ))
    expect_true("my_layout" %in% names(CellGraphs(se_layout)[[1]]@layout))
    expect_equal(dim(CellGraphs(se_layout)[[1]]@layout[["my_layout"]]), c(2470, 2))
  })


  test_that("ComputeLayout fails when invalid input is provided", {
    # Invalid object
    expect_error(ComputeLayout("Invalid"))
    e <- tryCatch(se[["mpxCells"]] %>% ComputeLayout(layout_method = "Invalid", verbose = FALSE), error = function(e) e, silent = TRUE)
    expect_s3_class(e, "simpleError")

    # Invalid combination of normalize_layout and project_on_unit_sphere
    expect_error(
      se[["mpxCells"]] %>% ComputeLayout(normalize_layout = TRUE, project_on_unit_sphere = TRUE, dim = 3)
    )

    # Invalid combination of dim and project_on_unit_sphere
    expect_error(
      se[["mpxCells"]] %>% ComputeLayout(dim = 2, project_on_unit_sphere = TRUE),
      "Projecting onto a unit sphere is only possible for 3D layouts"
    )

    custom_layout_fkn <- function(g) {
      return("Invalid")
    }

    cg <- CellGraphs(se)[[colnames(se)[1]]]

    # tbl_graph
    expect_error(layout <- ComputeLayout(cg@cellgraph, custom_layout_function = custom_layout_fkn))

    # CellGraph
    expect_error(cg_layout <- ComputeLayout(cg, custom_layout_function = custom_layout_fkn))

    # CellGraphAssay
    expect_error(cg_assay_layout <- ComputeLayout(se[["mpxCells"]], custom_layout_function = custom_layout_fkn))

    # Seurat
    expect_error(se_layout <- ComputeLayout(se, custom_layout_function = custom_layout_fkn))
  })
}
