options(Seurat.object.assay.version = "v3")

pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
cg <- CellGraphs(seur_obj)[[1]]

test_that("color_by_marker works as expected", {
  # Test with one marker
  expect_no_error(cg_colored <- color_by_marker(cg, markers = "CD9"))
  expect_true("val" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_true("color" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_equal(
    cg_colored@cellgraph %>% pull(val) %>% head(),
    c(
      `GGTATATATTTTAAGTAGTTTAGTA-A` = 0, `ATTATGTTGTTGTATGTTTATTATT-A` = 0,
      `TCCGTGGTTTGATACTCGGGAATTT-A` = 0, `AGTGTAAGAGGTTGTTTCTTAGAAA-A` = 0.167054084663166,
      `GAGCAGACAATGGCGCTTAGCTAAA-A` = 0, `CAATTTTGACCTAGTTGTGGCCAAG-A` = 0
    )
  )

  # Test with two markers
  expect_no_error(cg_colored <- color_by_marker(cg, markers = c("CD9", "ACTB")))
  expect_true("val" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_true("color" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_equal(
    cg_colored@cellgraph %>% pull(color) %>% head(),
    c("#440154", "#440154", "#440154", "#440154", "#440154", "#440154")
  )

  # Test with smooth counts
  expect_no_error(cg_colored <- color_by_marker(cg, markers = "CD9", smooth_counts = TRUE))
  expect_true("val" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_true("color" %in% names(cg_colored@cellgraph %>% as_tibble()))
  expect_equal(
    cg_colored@cellgraph %>% pull(color) %>% head(),
    c("#440154", "#440154", "#440154", "#3C4F8A", "#440154", "#440154")
  )

  # Test with nNodes
  expect_no_error(g_colored <- color_by_marker(cg, markers = "CD9", nNodes = 100))
  expect_equal(length(g_colored@cellgraph), 100)

  # Test with normalize FALSE
  expect_no_error(cg_colored <- color_by_marker(cg, markers = "CD9", normalize = FALSE))
  expect_equal(
    cg_colored@cellgraph %>% pull(color) %>% head(),
    c("#440154", "#440154", "#440154", "#35B779", "#440154", "#440154")
  )

  # Test with trim quantiles
  cg_colored <- color_by_marker(cg, markers = "CD9", trim_quantiles = c(0, 0.5))
  expect_equal(
    cg_colored@cellgraph %>% pull(color) %>% head(),
    c("#21908D", "#21908D", "#21908D", "#21908D", "#21908D", "#21908D")
  )

  # Test with mode = "sum"
  expect_no_error(cg_colored <- color_by_marker(cg, markers = c("CD9", "ACTB"), mode = "sum"))
  expect_equal(
    cg_colored@cellgraph %>% pull(color) %>% head(),
    c("#440154", "#440154", "#440154", "#3C4F8A", "#440154", "#440154")
  )
})


test_that("color_by_marker fails when invalid input is provided", {
  expect_error(cg_colored <- color_by_marker(cg, markers = "Invalid"), "Invalid missing from count matrix")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", smooth_counts = "Invalid"), "'smooth' must be TRUE or FALSE")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", mode = "Invalid"))
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", normalize = "Invalid"), "'normalize' must be TRUE or FALSE")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", trim_quantiles = c(-1, 1)), "'trim_quantiles' must be betweeen 0 and 1")
  expect_error(cg_colored <- color_by_marker(cg, markers = "CD9", nNodes = 1e6), "'nNodes' must be a positive value smaller than the number of nodes in the graph.")
})
