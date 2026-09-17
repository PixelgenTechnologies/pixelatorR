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

  expect_error(
    SeuratObject::FetchData(cg, vars = c("CD3", "not_a_variable")),
    "The following requested variables were not found"
  )
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

test_that("FetchData.CellGraph can add marker labels from one-hot counts", {
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
    counts = one_hot
  )

  fd <- SeuratObject::FetchData(cg_one_hot, vars = NULL, add_marker = TRUE)
  expect_equal(colnames(fd), "marker")
  expect_equal(fd$marker, c("CD3", "CD8", "CD4", "HLA-DR"))
  expect_equal(rownames(fd), node_names)

  fd_vars <- SeuratObject::FetchData(
    cg_one_hot,
    vars = "node_type",
    add_marker = TRUE,
    clean = FALSE
  )
  expect_equal(colnames(fd_vars), c("marker", "node_type"))
  expect_equal(fd_vars$marker, c("CD3", "CD8", "CD4", "HLA-DR"))

  expect_error(
    SeuratObject::FetchData(cg_one_hot, vars = "marker", add_marker = TRUE),
    "marker"
  )

  cg_empty <- CreateCellGraphObject(cellgraph = bipart_graph)
  fd_empty <- SeuratObject::FetchData(cg_empty, vars = NULL, add_marker = TRUE)
  expect_true(all(is.na(fd_empty$marker)))

  cg_na <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = one_hot,
    meta.data = data.frame(score = c(1, NA, 2, NA), row.names = node_names)
  )
  expect_warning(
    fd_clean <- SeuratObject::FetchData(
      cg_na,
      vars = "score",
      add_marker = TRUE,
      clean = TRUE
    ),
    "missing data for vars requested"
  )
  expect_equal(rownames(fd_clean), node_names[c(1, 3)])
  expect_equal(fd_clean$score, c(1, 2))
  expect_equal(fd_clean$marker, c("CD3", "CD4"))
})

test_that("FetchData.CellGraphList forwards add_marker", {
  one_hot <- Matrix::sparseMatrix(
    i = seq_len(4),
    j = c(1L, 1L, 2L, 2L),
    x = 1,
    dims = c(4, 4),
    dimnames = list(node_names, c("CD3", "CD4", "CD8", "HLA-DR"))
  )
  one_hot <- as(one_hot, "dgCMatrix")
  node_names_2 <- paste0("m", 1:4)
  bipart_graph_2 <- tidygraph::tbl_graph(
    nodes = data.frame(
      name = node_names_2,
      node_type = c("umi1", "umi1", "umi2", "umi2"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = c(1L, 2L, 3L), to = c(2L, 3L, 4L))
  )
  attr(bipart_graph_2, "type") <- "bipartite"
  one_hot_2 <- one_hot
  dimnames(one_hot_2) <- list(node_names_2, colnames(one_hot))
  cgl_one_hot <- CreateCellGraphList(list(
    cell_1 = CreateCellGraphObject(cellgraph = bipart_graph, counts = one_hot),
    cell_2 = CreateCellGraphObject(cellgraph = bipart_graph_2, counts = one_hot_2)
  ))

  fd <- SeuratObject::FetchData(cgl_one_hot, vars = "node_type", add_marker = TRUE)
  expect_equal(colnames(fd), c("component", "marker", "node_type"))
  expect_equal(fd$marker, rep(c("CD3", "CD3", "CD4", "CD4"), 2))
  expect_equal(rownames(fd), c(node_names, node_names_2))

  expect_error(
    SeuratObject::FetchData(cgl_one_hot, vars = "marker", add_marker = TRUE),
    "marker"
  )

  cgl_clean <- CreateCellGraphList(list(
    cell_1 = CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = one_hot,
      meta.data = data.frame(score = c(1, 2, 3, 4), row.names = node_names)
    ),
    cell_2 = CreateCellGraphObject(cellgraph = bipart_graph_2, counts = one_hot_2)
  ))
  expect_warning(
    expect_warning(
      fd_list_clean <- SeuratObject::FetchData(
        cgl_clean,
        vars = "score",
        add_marker = TRUE,
        clean = TRUE
      ),
      "missing from some graphs"
    ),
    "missing data for vars requested"
  )
  expect_equal(unique(fd_list_clean$component), "cell_1")
  expect_equal(nrow(fd_list_clean), 4)
  expect_true("marker" %in% colnames(fd_list_clean))
})

cg_no_cluster <- CreateCellGraphObject(
  cellgraph = bipart_graph,
  counts = counts
)
cgl <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg_no_cluster))

test_that("FetchData.CellGraphList works as expected", {
  expect_error(
    SeuratObject::FetchData(cgl, vars = c("CD3", "cluster", "node_type")),
    "duplicated"
  )

  fd_one <- SeuratObject::FetchData(cgl, vars = "CD3", cells = "cell_2")
  expect_s3_class(fd_one, "data.frame")
  expect_false(inherits(fd_one, "tbl_df"))
  expect_equal(unique(fd_one$component), "cell_2")
  expect_equal(nrow(fd_one), 4)
  expect_equal(rownames(fd_one), node_names)

  expect_error(
    SeuratObject::FetchData(cgl, vars = c("x", "CD3"), cells = "cell_1", clean = FALSE),
    "The following requested variables were not found"
  )

  expect_error(
    SeuratObject::FetchData(cgl, vars = c("CD3", "Invalid"), cells = "cell_1"),
    "The following requested variables were not found"
  )

  fd_hyphen <- SeuratObject::FetchData(cgl, vars = "HLA-DR", cells = "cell_1")
  expect_equal(colnames(fd_hyphen), c("component", "HLA-DR"))

  expect_error(SeuratObject::FetchData(cgl, vars = "component"), "cannot include")
})

test_that("FetchData.CellGraphList omits a missing layer on some graphs", {
  node_names_2 <- paste0("m", 1:4)
  bipart_graph_2 <- tidygraph::tbl_graph(
    nodes = data.frame(
      name = node_names_2,
      node_type = c("umi1", "umi1", "umi2", "umi2"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = c(1L, 2L, 3L), to = c(2L, 3L, 4L))
  )
  attr(bipart_graph_2, "type") <- "bipartite"
  counts_2 <- counts
  dimnames(counts_2) <- list(node_names_2, colnames(counts))
  cgl_lps <- CreateCellGraphList(list(
    cell_1 = CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      layers = list(lps = layer_mat)
    ),
    cell_2 = CreateCellGraphObject(cellgraph = bipart_graph_2, counts = counts_2)
  ))
  expect_warning(
    fd <- SeuratObject::FetchData(cgl_lps, vars = "f1", layer = "lps"),
    "missing from some graphs"
  )
  expect_equal(colnames(fd), c("component", "f1"))
  expect_equal(fd$f1[1:4], unname(layer_mat[, "f1"]))
  expect_true(all(is.na(fd$f1[5:8])))

  meta_2 <- data.frame(
    cluster = c("c", "c", "d", "d"),
    row.names = node_names_2,
    stringsAsFactors = FALSE
  )
  cgl_mixed <- CreateCellGraphList(list(
    cell_1 = CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      layers = list(lps = layer_mat),
      meta.data = meta
    ),
    cell_2 = CreateCellGraphObject(
      cellgraph = bipart_graph_2,
      counts = counts_2,
      meta.data = meta_2
    )
  ))
  expect_warning(
    fd_mixed <- SeuratObject::FetchData(
      cgl_mixed,
      vars = c("cluster", "f1"),
      layer = "lps"
    ),
    "missing from some graphs"
  )
  expect_equal(fd_mixed$cluster, c(meta$cluster, meta_2$cluster))
  expect_equal(fd_mixed$f1[1:4], unname(layer_mat[, "f1"]))
  expect_true(all(is.na(fd_mixed$f1[5:8])))
})

test_that("FetchData.CellGraphList keeps one row per node when no vars are requested", {
  node_names_2 <- paste0("m", 1:4)
  bipart_graph_2 <- tidygraph::tbl_graph(
    nodes = data.frame(
      name = node_names_2,
      node_type = c("umi1", "umi1", "umi2", "umi2"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = c(1L, 2L, 3L), to = c(2L, 3L, 4L))
  )
  attr(bipart_graph_2, "type") <- "bipartite"
  cgl_unique <- CreateCellGraphList(list(
    cell_1 = cg,
    cell_2 = CreateCellGraphObject(cellgraph = bipart_graph_2)
  ))

  expect_error(
    SeuratObject::FetchData(cgl_unique, vars = "Invalid"),
    "The following requested variables were not found"
  )

  fd_empty <- SeuratObject::FetchData(cgl_unique, vars = NULL)
  expect_equal(colnames(fd_empty), "component")
  expect_equal(nrow(fd_empty), 8)
})

test_that("FetchData.CellGraphList binds graphs with unique node IDs", {
  node_names_2 <- paste0("m", 1:4)
  bipart_graph_2 <- tidygraph::tbl_graph(
    nodes = data.frame(
      name = node_names_2,
      node_type = c("umi1", "umi1", "umi2", "umi2"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = c(1L, 2L, 3L), to = c(2L, 3L, 4L))
  )
  attr(bipart_graph_2, "type") <- "bipartite"
  counts_2 <- counts
  dimnames(counts_2) <- list(node_names_2, colnames(counts))
  meta_2 <- data.frame(
    cluster = c("c", "c", "d", "d"),
    row.names = node_names_2,
    stringsAsFactors = FALSE
  )
  cg_2 <- CreateCellGraphObject(
    cellgraph = bipart_graph_2,
    counts = counts_2,
    meta.data = meta_2
  )
  cgl_unique <- CreateCellGraphList(list(cell_1 = cg, cell_2 = cg_2))

  fd <- SeuratObject::FetchData(cgl_unique, vars = c("CD3", "cluster", "node_type"))
  expect_equal(colnames(fd), c("component", "CD3", "cluster", "node_type"))
  expect_equal(nrow(fd), 8)
  expect_equal(rownames(fd), c(node_names, node_names_2))
  expect_equal(unique(fd$component), c("cell_1", "cell_2"))
  expect_equal(fd$CD3, c(as.numeric(counts[, "CD3"]), as.numeric(counts_2[, "CD3"])))
  expect_equal(fd$cluster, c(meta$cluster, meta_2$cluster))
  expect_equal(fd$node_type, rep(c("umi1", "umi1", "umi2", "umi2"), 2))

  cgl_missing <- CreateCellGraphList(list(
    cell_1 = cg,
    cell_2 = CreateCellGraphObject(cellgraph = bipart_graph_2, counts = counts_2)
  ))
  expect_warning(
    expect_warning(
      fd_clean <- SeuratObject::FetchData(cgl_missing, vars = "cluster", clean = TRUE),
      "missing from some graphs"
    ),
    "missing data for vars requested"
  )
  expect_equal(unique(fd_clean$component), "cell_1")
  expect_equal(nrow(fd_clean), 4)
})

test_that("FetchData.CellGraphList keeps classed columns missing from a graph", {
  node_names_2 <- paste0("m", 1:4)
  bipart_graph_2 <- tidygraph::tbl_graph(
    nodes = data.frame(
      name = node_names_2,
      node_type = c("umi1", "umi1", "umi2", "umi2"),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = c(1L, 2L, 3L), to = c(2L, 3L, 4L))
  )
  attr(bipart_graph_2, "type") <- "bipartite"
  counts_2 <- counts
  dimnames(counts_2) <- list(node_names_2, colnames(counts))
  cg_factor <- CreateCellGraphObject(
    cellgraph = bipart_graph_2,
    counts = counts_2,
    meta.data = data.frame(
      grp = factor(c("a", "a", "b", "b"), levels = c("a", "b", "c")),
      day = as.Date("2020-01-01") + 0:3,
      row.names = node_names_2
    )
  )
  # The graph without the variables comes first, so the NA fill would
  # otherwise set the column type for the whole result
  cgl_classed <- CreateCellGraphList(list(cell_1 = cg_no_cluster, cell_2 = cg_factor))

  expect_warning(
    fd <- SeuratObject::FetchData(cgl_classed, vars = c("grp", "day"), clean = FALSE),
    "missing from some graphs"
  )
  expect_s3_class(fd$grp, "factor")
  expect_equal(levels(fd$grp), c("a", "b", "c"))
  expect_true(all(is.na(fd$grp[1:4])))
  expect_equal(as.character(fd$grp[5:8]), c("a", "a", "b", "b"))
  expect_s3_class(fd$day, "Date")
  expect_true(all(is.na(fd$day[1:4])))
  expect_equal(fd$day[5:8], as.Date("2020-01-01") + 0:3)

  cg_factor_a <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = counts,
    meta.data = data.frame(
      grp = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
      row.names = node_names
    )
  )
  cg_factor_c <- CreateCellGraphObject(
    cellgraph = bipart_graph_2,
    counts = counts_2,
    meta.data = data.frame(
      grp = factor(c("c", "c", "d", "d"), levels = c("c", "d")),
      row.names = node_names_2
    )
  )
  cgl_levels <- CreateCellGraphList(list(cell_1 = cg_factor_a, cell_2 = cg_factor_c))
  fd_levels <- SeuratObject::FetchData(cgl_levels, vars = "grp", clean = FALSE)
  expect_s3_class(fd_levels$grp, "factor")
  expect_equal(levels(fd_levels$grp), c("a", "b", "c", "d"))
  expect_equal(as.character(fd_levels$grp), c("a", "a", "b", "b", "c", "c", "d", "d"))
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
