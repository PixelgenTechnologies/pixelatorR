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

layer_mat <- matrix(
  seq_len(8),
  nrow = 4,
  ncol = 2,
  dimnames = list(node_names, c("f1", "f-2"))
)
meta <- data.frame(
  cluster = c("a", "a", "b", "b"),
  row.names = node_names,
  stringsAsFactors = FALSE
)
dr <- CreateNodeDimReducObject(
  embeddings = matrix(
    seq_len(8),
    nrow = 4,
    ncol = 2,
    dimnames = list(node_names, NULL)
  ),
  key = "PC_",
  method = "pca"
)

cg <- CreateCellGraphObject(
  cellgraph = bipart_graph,
  counts = counts,
  layers = list(data = layer_mat),
  meta.data = meta,
  reductions = list(pca = dr)
)

test_that("FetchData.CellGraph works as expected", {
  expect_no_error(fd <- SeuratObject::FetchData(cg, vars = c("CD3", "cluster", "PC_1", "node_type")))
  expect_s3_class(fd, "data.frame")
  expect_equal(nrow(fd), 4)
  expect_equal(rownames(fd), node_names)
  expect_equal(colnames(fd), c("CD3", "cluster", "PC_1", "node_type"))
  expect_equal(fd$CD3, as.numeric(counts[, "CD3"]))
  expect_equal(fd$cluster, meta$cluster)
  expect_equal(fd$PC_1, as.numeric(SeuratObject::Embeddings(cg, reduction = "pca")[, "PC_1"]))
  expect_equal(fd$node_type, c("umi1", "umi1", "umi2", "umi2"))

  fd_layer <- SeuratObject::FetchData(cg, vars = "f1", layer = "data")
  expect_equal(fd_layer$f1, unname(layer_mat[, "f1"]))

  expect_warning(fd_alt <- SeuratObject::FetchData(cg, vars = "f1"))
  expect_equal(fd_alt$f1, unname(layer_mat[, "f1"]))

  fd_cells <- SeuratObject::FetchData(cg, vars = "CD3", cells = node_names[1:2])
  expect_equal(rownames(fd_cells), node_names[1:2])
  expect_equal(fd_cells$CD3, as.numeric(counts[1:2, "CD3"]))

  fd_idx <- SeuratObject::FetchData(cg, vars = "CD3", cells = 1:2)
  expect_equal(rownames(fd_idx), node_names[1:2])

  expect_no_warning(
    fd_dup <- SeuratObject::FetchData(cg, vars = "CD3", cells = node_names[c(1, 1, 2)])
  )
  expect_equal(rownames(fd_dup), node_names[1:2])

  empty <- SeuratObject::FetchData(cg, vars = NULL)
  expect_equal(nrow(empty), 4)
  expect_equal(ncol(empty), 0)
})

test_that("FetchData.CellGraph keeps non-syntactic marker names", {
  expect_no_error(fd <- SeuratObject::FetchData(cg, vars = c("HLA-DR", "CD3")))
  expect_equal(colnames(fd), c("HLA-DR", "CD3"))
  expect_equal(fd[["HLA-DR"]], as.numeric(counts[, "HLA-DR"]))

  fd_layer <- SeuratObject::FetchData(cg, vars = "f-2", layer = "data")
  expect_equal(colnames(fd_layer), "f-2")
  expect_equal(fd_layer[["f-2"]], unname(layer_mat[, "f-2"]))

  fd_mixed <- SeuratObject::FetchData(cg, vars = c("HLA-DR", "cluster"))
  expect_equal(colnames(fd_mixed), c("HLA-DR", "cluster"))
  expect_equal(fd_mixed[["HLA-DR"]], as.numeric(counts[, "HLA-DR"]))
  expect_equal(fd_mixed$cluster, meta$cluster)
})

test_that("FetchData.CellGraph fails with invalid input", {
  expect_error(SeuratObject::FetchData(cg, vars = "not_a_variable"))
  expect_error(SeuratObject::FetchData(cg, vars = "CD3", cells = "missing_node"))
  expect_error(SeuratObject::FetchData(cg, vars = "CD3", layer = "missing_layer"))
  expect_warning(SeuratObject::FetchData(cg, vars = c("CD3", "missing")))
  expect_warning(
    SeuratObject::FetchData(cg, vars = "CD3", cells = c(node_names[1], "missing_node")),
    "1 node not present"
  )
})

cg_no_cluster <- CreateCellGraphObject(
  cellgraph = bipart_graph,
  counts = counts
)
cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg_no_cluster))

test_that("FetchData.CellGraphList works as expected", {
  expect_no_error(fd <- SeuratObject::FetchData(cgl, vars = c("CD3", "cluster", "node_type")))
  expect_s3_class(fd, "tbl_df")
  expect_equal(colnames(fd), c("component", "CD3", "cluster", "node_type"))
  expect_equal(nrow(fd), 8)
  expect_equal(unique(fd$component), c("cell_1", "cell_2"))
  expect_equal(fd$CD3, rep(as.numeric(counts[, "CD3"]), 2))
  expect_equal(fd$cluster[1:4], meta$cluster)
  expect_true(all(is.na(fd$cluster[5:8])))
  expect_equal(fd$node_type, rep(c("umi1", "umi1", "umi2", "umi2"), 2))

  fd_one <- SeuratObject::FetchData(cgl, vars = "CD3", cells = "cell_2")
  expect_equal(unique(fd_one$component), "cell_2")
  expect_equal(nrow(fd_one), 4)

  fd_xyz <- SeuratObject::FetchData(cgl, vars = c("x", "CD3"), clean = FALSE)
  expect_equal(colnames(fd_xyz), c("component", "x", "CD3"))
  expect_true(all(is.na(fd_xyz$x)))
  expect_equal(fd_xyz$CD3, rep(as.numeric(counts[, "CD3"]), 2))

  fd_hyphen <- SeuratObject::FetchData(cgl, vars = "HLA-DR")
  expect_equal(colnames(fd_hyphen), c("component", "HLA-DR"))

  expect_error(SeuratObject::FetchData(cgl, vars = "component"), "cannot include")

  fd_cluster <- SeuratObject::FetchData(cgl, vars = "cluster")
  expect_equal(nrow(fd_cluster), 8)
  expect_equal(unique(fd_cluster$component), c("cell_1", "cell_2"))
  expect_true(all(is.na(fd_cluster$cluster[5:8])))

  expect_warning(
    fd_clean <- SeuratObject::FetchData(cgl, vars = "cluster", clean = TRUE)
  )
  expect_equal(unique(fd_clean$component), "cell_1")
  expect_equal(nrow(fd_clean), 4)
})

test_that("FetchData.CellGraphList validates loaded CellGraphs", {
  cgl_mixed <- CreateCellGraphList(list(cell_1 = cg, cell_2 = NULL))
  expect_no_error(fd <- SeuratObject::FetchData(cgl_mixed, vars = "CD3"))
  expect_equal(unique(fd$component), "cell_1")
  expect_error(SeuratObject::FetchData(cgl_mixed, vars = "CD3", cells = c("cell_1", "cell_2")))
  expect_error(SeuratObject::FetchData(cgl_mixed, vars = "CD3", cells = "missing_cell"))

  cgl_empty <- CreateCellGraphList(list(cell_1 = NULL, cell_2 = NULL))
  expect_error(SeuratObject::FetchData(cgl_empty, vars = "CD3"))
})
