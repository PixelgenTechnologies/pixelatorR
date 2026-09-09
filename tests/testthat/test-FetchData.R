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
  1:12,
  nrow = 4,
  ncol = 3,
  sparse = TRUE,
  dimnames = list(node_names, c("CD3", "CD4", "CD8"))
)
counts <- as(counts, "dgCMatrix")

layer_mat <- matrix(
  seq_len(8),
  nrow = 4,
  ncol = 2,
  dimnames = list(node_names, c("f1", "f2"))
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
  expect_equal(fd_layer$f1, layer_mat[, "f1"])

  expect_warning(fd_alt <- SeuratObject::FetchData(cg, vars = "f1"))
  expect_equal(fd_alt$f1, layer_mat[, "f1"])

  fd_cells <- SeuratObject::FetchData(cg, vars = "CD3", cells = node_names[1:2])
  expect_equal(rownames(fd_cells), node_names[1:2])
  expect_equal(fd_cells$CD3, as.numeric(counts[1:2, "CD3"]))

  fd_idx <- SeuratObject::FetchData(cg, vars = "CD3", cells = 1:2)
  expect_equal(rownames(fd_idx), node_names[1:2])

  empty <- SeuratObject::FetchData(cg, vars = NULL)
  expect_equal(nrow(empty), 4)
  expect_equal(ncol(empty), 0)
})

test_that("FetchData.CellGraph fails with invalid input", {
  expect_error(SeuratObject::FetchData(cg, vars = "not_a_variable"))
  expect_error(SeuratObject::FetchData(cg, vars = "CD3", cells = "missing_node"))
  expect_error(SeuratObject::FetchData(cg, vars = "CD3", layer = "missing_layer"))
  expect_warning(SeuratObject::FetchData(cg, vars = c("CD3", "missing")))
})
