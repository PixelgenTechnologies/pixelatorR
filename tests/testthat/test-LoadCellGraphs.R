options(pixelatorR.arrow_outdir = tempdir())

# Load example data as a Seurat object
pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE, return_cellgraphassay = TRUE)
seur_obj_merged <- merge(seur_obj, seur_obj, add.cell.ids = c("Sample1", "Sample2"))
cg_assay <- seur_obj[["mpxCells"]]
cg_assay_merged <- merge(cg_assay, cg_assay)

test_that("LoadCellGraphs works for Seurat objects", {
  # Single data set
  expect_no_error({seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])})
  expect_s4_class(seur_obj, "Seurat")

  # Merged data set
  seur_obj_merged <- LoadCellGraphs(seur_obj_merged, cells = colnames(seur_obj_merged)[1])
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

test_that("LoadCellGraphs works for FileSystemDataset", {
  expect_no_error(el <- ReadMPX_arrow_edgelist(pxl_file))

  # Single data set
  expect_no_error({g_list <- LoadCellGraphs(el, cells = colnames(seur_obj)[1:2])})
  expect_type(g_list, "list")
  expect_s4_class(g_list[[1]], "CellGraph")
  expect_equal(length(g_list), 2)
})
