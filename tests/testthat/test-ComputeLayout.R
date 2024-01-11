se <- ReadMPX_Seurat(system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR"),
                     overwrite = TRUE, return_cellgraphassay = TRUE)

se <- LoadCellGraphs(se, cells = colnames(se)[1:2])
se[["mpxCells"]] <- KeepLargestComponent(se[["mpxCells"]])

test_that("ComputeLayout works as expected", {

  # Seurat object
  expect_no_error(se <- se %>% ComputeLayout())

  # CellGraphAssay
  layout_method <- "pmds"
  expect_no_error(cg_assay <- se[["mpxCells"]] %>% ComputeLayout(layout_method = layout_method))
  expect_true(layout_method %in% names(cg_assay@cellgraphs$RCVCMP0000000@layout))
  expect_equal(c("x", "y"), colnames(cg_assay@cellgraphs$RCVCMP0000000@layout[[layout_method]]))

  # Test with three dimensions
  expect_no_error(cg_assay <- se[["mpxCells"]] %>% ComputeLayout(dim = 3))
  expect_equal(c("x", "y", "z"), colnames(cg_assay@cellgraphs$RCVCMP0000000@layout[[layout_method]]))

  # Test with normalize_layout
  expect_no_error(cg_assay <- se[["mpxCells"]] %>% ComputeLayout(dim = 3, layout_method = "pmds", normalize_layout = TRUE))
  max_radius <- cg_assay@cellgraphs$RCVCMP0000000@layout[[layout_method]] %>%
    mutate(across(x:z, ~.x^2)) %>%
    rowSums() %>%
    sqrt() %>%
    max()
  expect_equal(max_radius, 1)

  # Test with project_on_unit_sphere
  expect_no_error(cg_assay <- se[["mpxCells"]] %>% ComputeLayout(dim = 3, layout_method = "pmds", project_on_unit_sphere = TRUE))
  radii <- cg_assay@cellgraphs$RCVCMP0000000@layout[[layout_method]] %>%
    mutate(across(x:z, ~.x^2)) %>%
    rowSums() %>%
    sqrt() %>%
    max()
  expect_equal(unique(radii) %>% length(), 1)
})

test_that("ComputeLayout works as expected with a custom layout function", {

  custom_layout_fkn <- graphlayouts::layout_with_pmds

  cg <- CellGraphs(se)[[colnames(se)[1]]]

  # tbl_graph
  expect_no_error(layout <- ComputeLayout(cg@cellgraph, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
  expect_equal(dim(layout), c(3525, 2))

  # CellGraph
  expect_no_error(cg_layout <- ComputeLayout(cg, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
  expect_equal(names(cg_layout@layout), "custom")
  expect_equal(dim(cg_layout@layout[["custom"]]), c(3525, 2))

  # CellGraphAssay
  expect_no_error(cg_assay_layout <- ComputeLayout(se[["mpxCells"]], custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
  expect_equal(names(CellGraphs(cg_assay_layout)[[1]]@layout), "custom")
  expect_equal(dim(CellGraphs(cg_assay_layout)[[1]]@layout[["custom"]]), c(3525, 2))

  # Seurat
  expect_no_error(se_layout <- ComputeLayout(se, custom_layout_function = custom_layout_fkn, custom_layout_function_args = list(pivots = 100)))
  expect_equal(names(CellGraphs(se_layout)[[1]]@layout), "custom")
  expect_equal(dim(CellGraphs(se_layout)[[1]]@layout[["custom"]]), c(3525, 2))

  # Test with new layout name
  expect_no_error(se_layout <- ComputeLayout(se, custom_layout_function = custom_layout_fkn,
                                             custom_layout_function_args = list(pivots = 100), custom_layout_name = "my_layout"))
  expect_equal(names(CellGraphs(se_layout)[[1]]@layout), "my_layout")
  expect_equal(dim(CellGraphs(se_layout)[[1]]@layout[["my_layout"]]), c(3525, 2))

})

test_that("ComputeLayout fails with an invalid input", {

  custom_layout_fkn <- function(g) return("Invalid")

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

test_that("ComputeLayout fails when invalid input is provided", {
  expect_error(ComputeLayout("Invalid"))
  e <- tryCatch(se[["mpxCells"]] %>% ComputeLayout(layout_method = "Invalid", verbose = FALSE), error = function(e) e, silent = TRUE)
  expect_s3_class(e, "simpleError")
})
