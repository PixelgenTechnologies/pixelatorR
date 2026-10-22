mpx_pxl_file <- minimal_mpx_pxl_file()
pna_pxl_file <- minimal_pna_pxl_file()

# Build a minimal MPX PXL file with precomputed layouts
pxl_file_new <- fs::path_temp("test.pxl")
seur_obj_mpx_subset <- ReadMPX_Seurat(mpx_pxl_file, verbose = FALSE)
seur_obj_mpx_subset <- subset(seur_obj_mpx_subset, cells = colnames(seur_obj_mpx_subset)[1:2]) |>
  LoadCellGraphs(verbose = FALSE) |>
  ComputeLayout(dim = 3, layout_method = "pmds", verbose = FALSE)
WriteMPX_pxl_file(seur_obj_mpx_subset, pxl_file_new, export_layouts = TRUE, overwrite = TRUE)

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expected_class_mpx <- ifelse(assay_version == "v3", "CellGraphAssay", "CellGraphAssay5")
  expected_class_pna <- ifelse(assay_version == "v3", "PNAAssay", "PNAAssay5")

  seur_obj_mpx <- ReadMPX_Seurat(mpx_pxl_file, verbose = FALSE)
  seur_obj_mpx_precomputed <- ReadMPX_Seurat(pxl_file_new, verbose = FALSE)
  cg_assay <- seur_obj_mpx[["mpxCells"]]
  cg_assay_split <- lapply(list(c(1, 2), c(3, 4, 5)), function(inds) {
    suppressWarnings(subset(cg_assay, cells = colnames(cg_assay)[inds]))
  })
  cg_assay_merged <- merge(
    cg_assay_split[[1]],
    cg_assay_split[-1],
    add.cell.ids = paste0("Sample", 1:2)
  )

  seur_obj_pna <- ReadPNA_Seurat(pna_pxl_file, verbose = FALSE)
  pna_assay <- seur_obj_pna[["PNA"]]
  pna_assay_split <- lapply(list(c(1, 2), c(3, 4, 5)), function(inds) {
    suppressWarnings(subset(pna_assay, cells = colnames(pna_assay)[inds]))
  })
  pna_assay_merged <- merge(
    pna_assay_split[[1]],
    pna_assay_split[-1],
    add.cell.ids = paste0("Sample", 1:2)
  )

  test_that("LoadCellGraphs works for Seurat objects", {
    expect_no_error({
      seur_obj_mpx <- LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], verbose = FALSE)
    })
    expect_s4_class(seur_obj_mpx, "Seurat")
    expect_no_error({
      seur_obj_pna <- LoadCellGraphs(seur_obj_pna, cells = colnames(seur_obj_pna)[1], verbose = FALSE)
    })
    expect_s4_class(seur_obj_pna, "Seurat")

    # Data set with pre-computed layouts (add_layouts is deprecated but still functional)
    lifecycle::expect_deprecated(
      seur_obj_mpx_precomputed <- LoadCellGraphs(
        seur_obj_mpx_precomputed,
        cells = colnames(seur_obj_mpx_precomputed)[1],
        add_layouts = TRUE,
        verbose = FALSE
      )
    )
    layouts <- seur_obj_mpx_precomputed[["mpxCells"]]@cellgraphs[[1]]@layout
    expect_equal(dim(layouts[[1]]), c(2470, 3))
    expect_equal(
      CellGraphData(seur_obj_mpx_precomputed[["mpxCells"]]@cellgraphs[[1]], slot = "nodes"),
      seur_obj_mpx_precomputed[["mpxCells"]]@cellgraphs[[1]]@cellgraph |> pull(name)
    )
    expect_lt(.row_names_info(layouts[[1]]), 0L)
    expect_equal(
      layouts[[1]] |> head() |> tibble::as_tibble(),
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

    lifecycle::expect_deprecated(
      seur_obj_pna <- LoadCellGraphs(
        seur_obj_pna,
        cells = colnames(seur_obj_pna)[1],
        add_layouts = TRUE,
        verbose = FALSE,
        force = TRUE
      )
    )
    layouts <- seur_obj_pna[["PNA"]]@cellgraphs[[1]]@layout
    expect_equal(dim(layouts[[1]]), c(43543, 3))
    expect_equal(
      CellGraphData(seur_obj_pna[["PNA"]]@cellgraphs[[1]], slot = "nodes"),
      seur_obj_pna[["PNA"]]@cellgraphs[[1]]@cellgraph |> pull(name)
    )
    expect_lt(.row_names_info(layouts[[1]]), 0L)
    expect_equal(
      layouts[[1]] |> head() |> tibble::as_tibble(),
      structure(
        list(
          x = c(
            -21845956.4386089,
            62412247.6490569,
            -51484652.1164919,
            40315293.3869398,
            18869639.379379,
            123838179.794506
          ),
          y = c(
            41934976.433566,
            -4536772.91578601,
            25394050.8184834,
            35106200.3000139,
            -27756606.072645,
            5526555.75287506
          ),
          z = c(
            20697597.0243981,
            23919096.3920098,
            -22291639.7321406,
            43829546.4013817,
            15284198.62568,
            -8333993.26120729
          )
        ),
        row.names = c(NA, -6L),
        class = c("tbl_df", "tbl", "data.frame")
      )
    )
  })

  test_that("LoadCellGraphs works for CellGraphAssay/PNAAssay objects", {
    # MPX projections
    for (proj in c("bipartite", "Anode", "linegraph")) {
      expect_no_error(
        cg_assay <- LoadCellGraphs(cg_assay, cells = colnames(cg_assay)[1], load_as = proj, force = TRUE, verbose = FALSE)
      )
      expect_s4_class(cg_assay, expected_class_mpx)
      expect_s4_class(CellGraphs(cg_assay)[[1]], "CellGraph")
      expect_equal(attr(CellGraphs(cg_assay)[[1]]@cellgraph, "type"), proj)
    }

    # PNA
    expect_no_error(pna_assay <- LoadCellGraphs(pna_assay, cells = colnames(pna_assay)[1], verbose = FALSE))
    expect_s4_class(pna_assay, expected_class_pna)
    expect_s4_class(CellGraphs(pna_assay)[[1]], "CellGraph")
    expect_equal(attr(CellGraphs(pna_assay)[[1]]@cellgraph, "type"), "bipartite")
  })

  test_that("LoadCellGraphs works for FileSystemDataset and merged assays", {
    expect_no_error(el <- ReadMPX_arrow_edgelist(mpx_pxl_file, verbose = FALSE))

    # Load 2 cells across 2 merged samples
    target_cells <- colnames(seur_obj_mpx)[c(1, 3)]
    expect_no_error(g_list <- LoadCellGraphs(el, cells = target_cells, data_type = "MPX", verbose = FALSE))
    expect_type(g_list, "list")
    expect_s4_class(g_list[[1]], "CellGraph")
    expect_equal(length(g_list), 2)

    merged_target_cells <- colnames(cg_assay_merged)[c(1, 3)]
    expect_no_error({
      g_list_merged <- LoadCellGraphs(cg_assay_merged, cells = merged_target_cells, verbose = FALSE)
    })
    expect_equal(names(CellGraphs(g_list_merged)[merged_target_cells]), merged_target_cells)

    for (i in seq_along(target_cells)) {
      cg_direct <- CellGraphData(g_list[[i]], slot = "cellgraph")
      cg_from_merged <- CellGraphData(CellGraphs(g_list_merged)[[merged_target_cells[i]]], slot = "cellgraph")
      expect_equal(igraph::gsize(cg_direct), igraph::gsize(cg_from_merged))
      expect_equal(length(cg_direct), length(cg_from_merged))
      expect_equal(attr(cg_direct, "type"), attr(cg_from_merged, "type"))
    }

    # Incremental loading (force = FALSE) and multi-sample merged PNA
    pna_merged_cells <- colnames(pna_assay_merged)
    pna_partial <- LoadCellGraphs(
      pna_assay_merged,
      cells = pna_merged_cells[c(1, 3)],
      verbose = FALSE
    )
    expect_s4_class(CellGraphs(pna_partial)[[pna_merged_cells[1]]], "CellGraph")
    expect_null(CellGraphs(pna_partial)[[pna_merged_cells[2]]])
    expect_s4_class(CellGraphs(pna_partial)[[pna_merged_cells[3]]], "CellGraph")

    attr(CellGraphs(pna_partial)[[pna_merged_cells[1]]]@cellgraph, "sentinel") <- "kept"
    expect_message(
      pna_incremental <- LoadCellGraphs(
        pna_partial,
        cells = pna_merged_cells[1:4],
        force = FALSE,
        verbose = TRUE
      ),
      "2 CellGraphs loaded. Loading remaining 2 CellGraphs."
    )
    expect_equal(
      attr(CellGraphs(pna_incremental)[[pna_merged_cells[1]]]@cellgraph, "sentinel"),
      "kept"
    )
    expect_s4_class(CellGraphs(pna_incremental)[[pna_merged_cells[2]]], "CellGraph")
    expect_s4_class(CellGraphs(pna_incremental)[[pna_merged_cells[4]]], "CellGraph")
    expect_null(CellGraphs(pna_incremental)[[pna_merged_cells[5]]])
  })

  test_that("LoadCellGraphs fails with invalid input", {
    expect_error(LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], load_as = "invalid"))
    expect_error(LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], assay = "invalid"))
    expect_error(LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], add_marker_counts = "invalid"))
    expect_error(suppressWarnings(LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], add_layouts = "invalid")))
    expect_error(LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], force = "invalid"))
    expect_error(LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], chunk_size = "invalid"))

    # LoadCellGraphs should fail with duplicated cell ids
    expect_error(LoadCellGraphs(seur_obj_mpx, cells = rep(colnames(seur_obj_mpx)[1], 2)))
  })
}

test_that("LoadCellGraphs warns when the deprecated add_layouts argument is used", {
  # Supplying add_layouts (with any value) triggers a deprecation warning ...
  lifecycle::expect_deprecated(
    LoadCellGraphs(seur_obj_pna, cells = colnames(seur_obj_pna)[1], add_layouts = TRUE, force = TRUE, verbose = FALSE),
    regexp = "add_layouts"
  )
  lifecycle::expect_deprecated(
    LoadCellGraphs(seur_obj_pna[["PNA"]], cells = colnames(seur_obj_pna)[1], add_layouts = FALSE, force = TRUE, verbose = FALSE),
    regexp = "add_layouts"
  )

  # ... but the deprecated argument still works
  lifecycle::expect_deprecated(
    seur_obj_pna <- LoadCellGraphs(seur_obj_pna, cells = colnames(seur_obj_pna)[1], add_layouts = TRUE, force = TRUE, verbose = FALSE)
  )
  expect_named(CellGraphs(seur_obj_pna)[[1]]@layout, "wpmds_3d")

  # No warning when add_layouts is not supplied
  expect_no_condition(
    LoadCellGraphs(seur_obj_pna, cells = colnames(seur_obj_pna)[1], force = TRUE, verbose = FALSE),
    class = "lifecycle_warning_deprecated"
  )
})

test_that(".attach_precomputed_layouts keeps distinct A/B coordinates", {
  suffixed <- c("umi1-A", "umi1-B", "umi2-A", "umi2-B")
  g <- tidygraph::tbl_graph(
    nodes = data.frame(
      name = suffixed,
      node_type = c("A", "B", "A", "B"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = c(1L, 3L), to = c(2L, 4L))
  )
  attr(g, "type") <- "bipartite"
  cg <- CreateCellGraphObject(cellgraph = g)

  # Native MPX: one row per UMI, A and B share coordinates
  unique_umi <- data.frame(
    name = c("umi2", "umi1"),
    x = c(20, 10),
    y = c(2, 1)
  )
  cg_unique <- pixelatorR:::.attach_precomputed_layouts(cg, list(pmds = unique_umi))
  expect_equal(cg_unique@layout$pmds$x, c(10, 10, 20, 20))
  expect_equal(cg_unique@layout$pmds$y, c(1, 1, 2, 2))

  # WriteMPX: one row per graph node with stripped names, distinct A/B coords
  written <- data.frame(
    name = c("umi1", "umi1", "umi2", "umi2"),
    x = c(10, 11, 20, 21),
    y = c(1, 1.1, 2, 2.1)
  )
  cg_written <- pixelatorR:::.attach_precomputed_layouts(cg, list(pmds = written))
  expect_equal(cg_written@layout$pmds$x, c(10, 11, 20, 21))
  expect_equal(cg_written@layout$pmds$y, c(1, 1.1, 2, 2.1))
})
