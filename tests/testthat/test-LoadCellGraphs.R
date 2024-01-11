# Load example data as a Seurat object
pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE, return_cellgraphassay = TRUE)
cg_assay <- seur_obj[["mpxCells"]]

test_that("LoadCellGraphs works for Seurat objects", {
  expect_no_error({seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])})
  expect_s4_class(seur_obj, "Seurat")
})

test_that("LoadCellGraphs works for CellGraphAssay objects", {

  # Load bipartite graph (default)
  expect_no_error({cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], force = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_s4_class(cg_assay@cellgraphs$RCVCMP0000000, "CellGraph")
  expect_equal(attr(cg_assay@cellgraphs$RCVCMP0000000@cellgraph, "type"), "bipartite")

  # Load A-node projected graph
  expect_no_error({cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = "Anode", force = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_s4_class(cg_assay@cellgraphs$RCVCMP0000000, "CellGraph")
  expect_equal(attr(cg_assay@cellgraphs$RCVCMP0000000@cellgraph, "type"), "Anode")

  # Load line graph
  expect_no_error({cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = "linegraph", force = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_s4_class(cg_assay@cellgraphs$RCVCMP0000000, "CellGraph")
  expect_equal(attr(cg_assay@cellgraphs$RCVCMP0000000@cellgraph, "type"), "linegraph")
})

