pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file)
pxl_file_new <- fs::path_temp("test.pxl")

# Load cell graphs and compute layouts
seur_obj <- seur_obj %>%
  LoadCellGraphs(cells = colnames(.)[1:2]) %>%
  ComputeLayout(dim = 3, layout_method = "pmds")

# Export layouts to a pxl file
seur_obj_subset <- subset(seur_obj, cells = colnames(seur_obj)[1:2])
WriteMPX_pxl_file(seur_obj_subset, pxl_file_new, export_layouts = TRUE, overwrite = TRUE)

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expected_class <- ifelse(assay_version == "v3", "CellGraphAssay", "CellGraphAssay5")

  # Load example data as a Seurat object
  pxl_file <- system.file("extdata/five_cells",
    "five_cells.pxl",
    package = "pixelatorR"
  )
  seur_obj <- ReadMPX_Seurat(pxl_file)
  seur_obj_merged <- merge(seur_obj, seur_obj, add.cell.ids = c("Sample1", "Sample2"))
  cg_assay <- seur_obj[["mpxCells"]]
  cg_assay_merged <- merge(cg_assay, cg_assay, add.cell.ids = c("A", "B"))

  seur_obj_precomputed <- ReadMPX_Seurat(pxl_file_new)

  cg_assay_split <- lapply(list(c(1, 2), c(3, 4, 5)), function(inds) {
    suppressWarnings({
      subset(cg_assay, cells = colnames(cg_assay)[inds])
    })
  })
  cg_assay_merged_big <- merge(cg_assay_split[[1]], cg_assay_split[-1], add.cell.ids = paste0("Sample", 1:2))

  test_that("LoadCellGraphs works for Seurat objects", {
    # Single data set
    expect_no_error({
      seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
    })
    expect_s4_class(seur_obj, "Seurat")

    # Data set with pre-computed layouts
    seur_obj_precomputed <- LoadCellGraphs(seur_obj_precomputed, cells = colnames(seur_obj)[1], load_layouts = TRUE)
    layouts <- seur_obj_precomputed[["mpxCells"]]@cellgraphs[[1]]@layout
    expect_equal(dim(layouts[[1]]), c(2470, 3))
    expect_equal(
      layouts[[1]] %>% head(),
      structure(
        list(
          x = c(
            -58.5281066499656,
            -47.4873876672378,
            -26.208082178209,
            60.4958073447489,
            -47.7818086804891,
            -55.4854489842481
          ),
          y = c(
            -4.41079380126333, -93.6673611597246,
            -5.05638704235889,
            1.8465975759391,
            44.5521164154744, -30.8407062275168
          ),
          z = c(
            -23.596714663682,
            5.60845902338665,
            36.2101781059961,
            56.9861893594522,
            1.57940672261431,
            -0.384017528242823
          )
        ),
        row.names = c(NA, -6L),
        class = c("tbl_df", "tbl", "data.frame")
      )
    )

    # Merged data set
    seur_obj_merged <- LoadCellGraphs(seur_obj_merged, cells = colnames(seur_obj_merged)[1])
  })

  test_that("LoadCellGraphs works for CellGraphAssay objects", {
    # Load bipartite graph (default)
    expect_no_error({
      cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], force = TRUE)
    })
    expect_s4_class(cg_assay, expected_class)
    expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
    expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), "bipartite")

    # Load A-node projected graph
    expect_no_error({
      cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = "Anode", force = TRUE)
    })
    expect_s4_class(cg_assay, expected_class)
    expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
    expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), "Anode")

    # Load line graph
    expect_no_error({
      cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = "linegraph", force = TRUE)
    })
    expect_s4_class(cg_assay, expected_class)
    expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
    expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), "linegraph")
  })

  test_that("LoadCellGraphs works for FileSystemDataset", {
    expect_no_error(el <- ReadMPX_arrow_edgelist(pxl_file))

    # Single data set
    expect_no_error({
      g_list <- LoadCellGraphs(el, cells = colnames(seur_obj))
    })
    expect_type(g_list, "list")
    expect_s4_class(g_list[[1]], "CellGraph")
    expect_equal(length(g_list), 5)

    # Big merged data
    expect_no_error({
      g_list_merged_big <- LoadCellGraphs(cg_assay_merged_big, cells = colnames(cg_assay_merged_big))
    })
    expect_equal(names(CellGraphs(g_list_merged_big)), colnames(cg_assay_merged_big))

    # Check that contents of g_list and g_list_merged_big are the same
    for (i in 1:2) {
      expect_equal(
        g_list[[i]] %>% CellGraphData(slot = "cellgraph") %>% igraph::gsize(),
        CellGraphs(g_list_merged_big)[[i]] %>% CellGraphData(slot = "cellgraph") %>% igraph::gsize()
      )
      expect_equal(
        g_list[[i]] %>% CellGraphData(slot = "cellgraph") %>% length(),
        CellGraphs(g_list_merged_big)[[i]] %>% CellGraphData(slot = "cellgraph") %>% length()
      )
      expect_equal(
        g_list[[i]] %>% CellGraphData(slot = "cellgraph") %>% attr(which = "type"),
        CellGraphs(g_list_merged_big)[[i]] %>% CellGraphData(slot = "cellgraph") %>% attr(which = "type")
      )
    }
  })
}
