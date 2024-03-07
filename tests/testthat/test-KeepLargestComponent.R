pxl_file <- system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR")
seur <- ReadMPX_Seurat(filename = pxl_file, return_cellgraphassay = TRUE, overwrite = TRUE)
seur <- LoadCellGraphs(seur, cells = colnames(seur)[1])

# Breaks graph
cg <- CellGraphs(seur)[[1]]
set.seed(123)
cg@cellgraph <- cg@cellgraph %E>%
  filter(from %in% sample(from, igraph::gsize(.) - 500, replace = FALSE))
seur@assays$mpxCells@cellgraphs[[1]] <- cg
cg_assay <- seur[["mpxCells"]]

test_that("KeepLargestComponent works as expected", {

  # Seurat
  expect_no_error(seur_large_component <- KeepLargestComponent(seur))
  expect_s4_class(seur_large_component, "Seurat")
  CellGraphs(seur) <- rep(list(NULL), ncol(seur)) %>% setNames(nm = colnames(seur))
  expect_no_error(seur <- KeepLargestComponent(seur))

  # CellGraphAssay
  expect_no_error(cg_assay_large_component <- KeepLargestComponent(cg_assay))
  expect_s4_class(cg_assay_large_component, "CellGraphAssay")
  cg_assay_large_component@cellgraphs <- rep(list(NULL), ncol(cg_assay)) %>% setNames(nm = colnames(cg_assay))
  expect_no_error(cg_assay_large_component <- KeepLargestComponent(cg_assay_large_component))

  # CellGraph
  expect_no_error(cg_large_component <- KeepLargestComponent(CellGraphs(cg_assay)[[1]]))
  expect_s4_class(cg_large_component, "CellGraph")

  # tbl_graph
  expect_no_error(edge_list_large_component <- KeepLargestComponent(CellGraphData(CellGraphs(cg_assay)[[1]], slot = "cellgraph"), verbose = FALSE))
  expect_s3_class(edge_list_large_component, "tbl_graph")
  expect_equal(igraph::gsize(edge_list_large_component), 5092)
  expect_equal(length(edge_list_large_component), 2429)
  expect_no_error(edge_list_large_component <- KeepLargestComponent(edge_list_large_component))
})

test_that("KeepLargestComponent fails when invalid input is provided", {
  expect_error(KeepLargestComponent(object = "Invalid"))
})
