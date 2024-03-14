options(pixelatorR.arrow_outdir = tempdir())

# Load example data as a Seurat object
pxl_file <- system.file("extdata/mock_data",
                        "mock_mpx_data.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE, return_cellgraphassay = TRUE)
seur_obj_merged <- merge(seur_obj, seur_obj, add.cell.ids = c("Sample1", "Sample2"))
cg_assay <- seur_obj[["mpxCells"]]
cg_assay_merged <- merge(cg_assay, cg_assay)


cg_assay_split <- lapply(list(c(1, 2), c(3, 4, 5)), function(inds) {
  subset(cg_assay, cells = colnames(cg_assay)[inds])
})
cg_assay_merged_big <- merge(cg_assay_split[[1]], cg_assay_split[-1], add.cell.ids = paste0("Sample", 1:2))

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
  expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
  expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), "bipartite")

  # Load A-node projected graph
  expect_no_error({cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = "Anode", force = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
  expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), "Anode")

  # Load line graph
  expect_no_error({cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = "linegraph", force = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
  expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), "linegraph")
})

test_that("LoadCellGraphs works for FileSystemDataset", {

  expect_no_error(el <- ReadMPX_arrow_edgelist(pxl_file))

  # Single data set
  expect_no_error({g_list <- LoadCellGraphs(el, cells = colnames(seur_obj))})
  expect_type(g_list, "list")
  expect_s4_class(g_list[[1]], "CellGraph")
  expect_equal(length(g_list), 5)

  # Big merged data
  expect_no_error({g_list_merged_big <- LoadCellGraphs(ArrowData(cg_assay_merged_big), cells = colnames(cg_assay_merged_big))})
  expect_equal(names(g_list_merged_big), colnames(cg_assay_merged_big))

  # Check data in arrow directory
  expect_equal(ArrowDir(cg_assay_merged_big) %>% list.files(),
               c("sample=Sample1", "sample=Sample2"))

  # Check that components in parquet files are correct
  for (i in 1:2) {
    expect_equal(
      ArrowData(cg_assay_merged_big) %>% filter(sample == paste0("Sample", i)) %>% collect() %>% pull(component) %>% unique() %>% as.character() %>% sort(),
      colnames(cg_assay_split[[i]]) %>% as.character() %>% sort()
    )
  }

  # Check that contents of g_list and g_list_merged_big are the same
  for (i in 1:2) {
    expect_equal(g_list[[i]] %>% CellGraphData(slot = "cellgraph") %>% igraph::gsize(),
                 g_list_merged_big[[i]] %>% CellGraphData(slot = "cellgraph") %>% igraph::gsize()
    )
    expect_equal(g_list[[i]] %>% CellGraphData(slot = "cellgraph") %>% length(),
                 g_list_merged_big[[i]] %>% CellGraphData(slot = "cellgraph") %>% length()
    )
    expect_equal(g_list[[i]] %>% CellGraphData(slot = "cellgraph") %>% attr(which = "type"),
                 g_list_merged_big[[i]] %>% CellGraphData(slot = "cellgraph") %>% attr(which = "type")
    )
  }
})
