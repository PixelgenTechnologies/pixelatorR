library(pixelatorR)
library(SeuratObject)

node_names <- paste0("n", 1:4)
bipart_graph <- tidygraph::tbl_graph(
  nodes = data.frame(
    name = node_names,
    node_type = c("umi1", "umi1", "umi2", "umi2"),
    stringsAsFactors = FALSE
  ),
  edges = data.frame(from = c(1L, 2L, 3L), to = c(2L, 3L, 4L))
)
attr(bipart_graph, "type") <- "bipartite"

counts <- Matrix::Matrix(
  1:16,
  nrow = 4,
  ncol = 4,
  sparse = TRUE,
  dimnames = list(node_names, c("CD3", "CD4", "CD8", "HLA-DR"))
)
counts <- as(counts, "dgCMatrix")

meta <- data.frame(
  cluster = factor(c("a", "a", "b", "b")),
  row.names = node_names,
  stringsAsFactors = FALSE
)

layout <- data.frame(
  x = c(0.1, 0.2, 0.3, 0.4),
  y = c(1.1, 1.2, 1.3, 1.4),
  z = c(2.1, 2.2, 2.3, 2.4),
  row.names = node_names
)

cg <- CreateCellGraphObject(
  cellgraph = bipart_graph,
  counts = counts,
  layout = list(wpmds_3d = layout),
  meta.data = meta
)

cg_no_cluster <- CreateCellGraphObject(
  cellgraph = bipart_graph,
  counts = counts,
  layout = list(wpmds_3d = layout)
)

test_that("FetchLayoutData.CellGraph works as expected", {
  expect_no_error(lyt <- FetchLayoutData(cg))
  expect_s3_class(lyt, "tbl_df")
  expect_equal(colnames(lyt), c("x", "y", "z"))
  expect_equal(nrow(lyt), 4)
  expect_equal(lyt$x, layout$x)
  expect_equal(lyt$y, layout$y)
  expect_equal(lyt$z, layout$z)

  expect_no_error(lyt_vars <- FetchLayoutData(cg, vars = c("CD3", "cluster", "node_type")))
  expect_equal(colnames(lyt_vars), c("x", "y", "z", "CD3", "cluster", "node_type"))
  expect_equal(lyt_vars$CD3, as.numeric(counts[, "CD3"]))
  expect_equal(lyt_vars$cluster, meta$cluster)
  expect_true(is.factor(lyt_vars$cluster))
  expect_equal(lyt_vars$node_type, c("umi1", "umi1", "umi2", "umi2"))
})

test_that("FetchLayoutData.CellGraph drops vars missing from the graph", {
  expect_warning(
    lyt <- FetchLayoutData(cg, vars = c("CD3", "missing_var")),
    "The following requested variables were not found"
  )
  expect_equal(colnames(lyt), c("x", "y", "z", "CD3"))
  expect_equal(lyt$CD3, as.numeric(counts[, "CD3"]))
  expect_false("missing_var" %in% colnames(lyt))

  expect_warning(
    lyt_all_missing <- FetchLayoutData(cg, vars = "not_a_variable"),
    "The following requested variables were not found"
  )
  expect_equal(colnames(lyt_all_missing), c("x", "y", "z"))
})

test_that("FetchLayoutData.CellGraph keeps non-syntactic marker names", {
  expect_no_error(lyt <- FetchLayoutData(cg, vars = "HLA-DR"))
  expect_equal(colnames(lyt), c("x", "y", "z", "HLA-DR"))
  expect_equal(lyt[["HLA-DR"]], as.numeric(counts[, "HLA-DR"]))
})

test_that("FetchLayoutData.CellGraph can add protein labels from one-hot counts", {
  one_hot <- Matrix::sparseMatrix(
    i = seq_len(4),
    j = c(1L, 3L, 2L, 4L),
    x = 1,
    dims = c(4, 4),
    dimnames = list(node_names, c("CD3", "CD4", "CD8", "HLA-DR"))
  )
  one_hot <- as(one_hot, "dgCMatrix")
  cg_one_hot <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = one_hot,
    layout = list(wpmds_3d = layout)
  )

  lyt <- FetchLayoutData(cg_one_hot, add_protein = TRUE)
  expect_equal(colnames(lyt), c("x", "y", "z", "protein"))
  expect_equal(lyt$protein, c("CD3", "CD8", "CD4", "HLA-DR"))

  lyt_vars <- FetchLayoutData(cg_one_hot, vars = "node_type", add_protein = TRUE)
  expect_equal(colnames(lyt_vars), c("x", "y", "z", "protein", "node_type"))

  cg_empty <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    layout = list(wpmds_3d = layout)
  )
  lyt_empty <- FetchLayoutData(cg_empty, add_protein = TRUE)
  expect_true(all(is.na(lyt_empty$protein)))
})

test_that("FetchLayoutData.CellGraphList forwards add_protein", {
  one_hot <- Matrix::sparseMatrix(
    i = seq_len(4),
    j = c(1L, 1L, 2L, 2L),
    x = 1,
    dims = c(4, 4),
    dimnames = list(node_names, c("CD3", "CD4", "CD8", "HLA-DR"))
  )
  one_hot <- as(one_hot, "dgCMatrix")
  cg_one_hot <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = one_hot,
    layout = list(wpmds_3d = layout)
  )
  cgl_one_hot <- CreateCellGraphList(list(cell_1 = cg_one_hot, cell_2 = cg_one_hot))
  lyt <- FetchLayoutData(cgl_one_hot, add_protein = TRUE)
  expect_equal(colnames(lyt), c("component", "x", "y", "z", "protein"))
  expect_equal(lyt$protein, rep(c("CD3", "CD3", "CD4", "CD4"), 2))
})

test_that("FetchLayoutData.CellGraph fails with invalid input", {
  expect_error(FetchLayoutData(cg, layout_method = "missing_layout"))
  expect_error(FetchLayoutData(cg, vars = "protein"))
  expect_error(FetchLayoutData(cg, vars = "x"))
  cg_no_layout <- CreateCellGraphObject(cellgraph = bipart_graph, counts = counts)
  expect_error(FetchLayoutData(cg_no_layout))
  cg_2d <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = counts,
    layout = list(wpmds = data.frame(x = 1:4, y = 1:4, row.names = node_names))
  )
  expect_error(FetchLayoutData(cg_2d, layout_method = "wpmds"))
})

cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg_no_cluster))

test_that("FetchLayoutData.CellGraphList works as expected", {
  expect_no_error(lyt <- FetchLayoutData(cgl, vars = c("CD3", "cluster")))
  expect_s3_class(lyt, "tbl_df")
  expect_equal(colnames(lyt), c("component", "x", "y", "z", "CD3", "cluster"))
  expect_equal(nrow(lyt), 8)
  expect_equal(unique(lyt$component), c("cell_1", "cell_2"))
  expect_equal(lyt$CD3, rep(as.numeric(counts[, "CD3"]), 2))
  expect_equal(lyt$cluster[1:4], meta$cluster)
  expect_true(all(is.na(lyt$cluster[5:8])))

  expect_no_error(lyt_one <- FetchLayoutData(cgl, cells = "cell_2", vars = "CD3"))
  expect_equal(unique(lyt_one$component), "cell_2")
  expect_equal(nrow(lyt_one), 4)
})

test_that("FetchLayoutData.CellGraphList drops vars missing from every graph", {
  expect_warning(
    lyt <- FetchLayoutData(cgl, vars = c("CD3", "missing_var")),
    "The following requested variables were not found"
  )
  expect_equal(colnames(lyt), c("component", "x", "y", "z", "CD3"))
  expect_false("missing_var" %in% colnames(lyt))
})

test_that("FetchLayoutData.CellGraphList validates loaded CellGraphs", {
  cgl_mixed <- CreateCellGraphList(list(cell_1 = cg, cell_2 = NULL))
  expect_no_error(lyt <- FetchLayoutData(cgl_mixed))
  expect_equal(unique(lyt$component), "cell_1")
  expect_error(FetchLayoutData(cgl_mixed, cells = c("cell_1", "cell_2")))
  expect_error(FetchLayoutData(cgl_mixed, cells = "missing_cell"))

  cgl_empty <- CreateCellGraphList(list(cell_1 = NULL, cell_2 = NULL))
  expect_error(FetchLayoutData(cgl_empty))
})

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
se <- LoadCellGraphs(se, cells = colnames(se)[1:2], add_layouts = TRUE, verbose = FALSE)
cells <- colnames(se)[1:2]

test_that("FetchLayoutData.PNAAssay and Seurat methods work as expected", {
  expect_no_error(lyt_assay <- FetchLayoutData(se[["PNA"]], cells = cells[1], vars = "B2M"))
  expect_s3_class(lyt_assay, "tbl_df")
  expect_true(all(c("component", "x", "y", "z", "B2M") %in% colnames(lyt_assay)))
  expect_equal(unique(lyt_assay$component), cells[1])
  expect_false(any(is.na(lyt_assay$x)))
  expect_false(any(is.na(lyt_assay$B2M)))

  expect_no_error(lyt_protein <- FetchLayoutData(se[["PNA"]], cells = cells[1], add_protein = TRUE))
  expect_true("protein" %in% colnames(lyt_protein))
  expect_equal(nrow(lyt_protein), nrow(lyt_assay))
  expect_true(is.character(lyt_protein$protein))
  expect_gt(length(unique(lyt_protein$protein)), 1)

  expect_no_error(lyt_seurat <- FetchLayoutData(se, cells = cells, vars = "B2M"))
  expect_equal(sort(unique(lyt_seurat$component)), sort(cells))
  expect_equal(nrow(lyt_seurat), sum(vapply(CellGraphs(se)[cells], function(x) length(CellGraphData(x, slot = "nodes")), integer(1))))

  expect_warning(
    lyt_invalid <- FetchLayoutData(se, cells = cells, vars = c("B2M", "Invalid"), add_protein = TRUE),
    "The following requested variables were not found"
  )
  expect_true(all(c("component", "x", "y", "z", "protein", "B2M") %in% colnames(lyt_invalid)))
  expect_false("Invalid" %in% colnames(lyt_invalid))
})

test_that("FetchLayoutData.PNAAssay and Seurat methods validate loaded CellGraphs", {
  se_unloaded <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
  expect_error(FetchLayoutData(se_unloaded))
  expect_error(FetchLayoutData(se_unloaded[["PNA"]]))
  expect_error(FetchLayoutData(se, cells = colnames(se)[3]))
})
