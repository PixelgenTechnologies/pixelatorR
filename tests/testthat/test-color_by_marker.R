pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
cg <- seur_obj[["mpxCells"]]@cellgraphs[[colnames(seur_obj)[1]]]

test_that("color_by_marker works as expected", {

  # Test with one marker
  expect_no_error(cg_colored <- color_by_marker(cg, markers = "CD9"))
  expect_true("val" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_true("color" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_equal(cg_colored@cellgraph %>% pull(val) %>% head() %>% unname(),
               c(0, 0, 0, 0, 0.09531018, 0.06453852))

  # Test with two markers
  expect_no_error(cg_colored <- color_by_marker(cg, markers = c("CD9", "ACTB")))
  expect_true("val" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_true("color" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_equal(cg_colored@cellgraph %>% pull(color) %>% head() %>% unname(),
               c("#440154", "#440154", "#440154", "#440154", "#440154", "#32658E"))

  # Test with smooth counts
  expect_no_error(cg_colored <- color_by_marker(cg, markers = "CD9", smooth_counts = TRUE))
  expect_true("val" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_true("color" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_equal(cg_colored@cellgraph %>% pull(color) %>% head() %>% unname(),
               c("#440154", "#440154", "#440154", "#440154", "#482071", "#481769"))

  # Test with nNodes
  expect_no_error(g_colored <- color_by_marker(cg, markers = "CD9", nNodes = 100))
  expect_equal(length(g_colored@cellgraph), 100)

  # Test with normalize FALSE
  expect_no_error(cg_colored <- color_by_marker(cg, markers = "CD9", normalize = FALSE))
  expect_equal(cg_colored@cellgraph %>% pull(color) %>% head() %>% unname(),
               c("#440154", "#440154", "#440154", "#440154", "#21908D", "#21908D"))

  # Test with trim quantiles
  cg_colored <- color_by_marker(cg, markers = "CD9", trim_quantiles = c(0, 0.5))
  expect_equal(cg_colored@cellgraph %>% pull(color) %>% head() %>% unname(),
               c("#21908D", "#21908D", "#21908D", "#21908D", "#21908D", "#21908D"))

  # Test with mode = "sum"
  expect_no_error(cg_colored <- color_by_marker(cg, markers = c("CD9", "ACTB"), mode = "sum"))
  expect_equal(cg_colored@cellgraph %>% pull(color) %>% head() %>% unname(),
               c("#440154", "#481466", "#440154", "#440154", "#482071", "#482979"))
})


test_that("color_by_marker fails when invalid input is provided", {

  expect_error(cg_colored <- color_by_marker(cg, markers = "Invalid"), "Invalid missing from count matrix")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", smooth_counts = "Invalid"), "'smooth' must be TRUE or FALSE")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", mode = "Invalid"))
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", normalize = "Invalid"), "'normalize' must be TRUE or FALSE")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", trim_quantiles = c(-1, 1)), "'trim_quantiles' must be betweeen 0 and 1")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", nNodes = 1e6), "'nNodes' must be a positive value smaller than the number of nodes in the graph.")
})
