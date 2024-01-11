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

test_that("ComputeLayout fails when invalid input is provided", {
  expect_error(ComputeLayout("Invalid"))
  e <- tryCatch(se[["mpxCells"]] %>% ComputeLayout(layout_method = "Invalid", verbose = FALSE), error = function(e) e, silent = TRUE)
  expect_s3_class(e, "simpleError")
})
