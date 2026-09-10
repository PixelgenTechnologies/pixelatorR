edge_list <-
  ReadMPX_item(
    system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
    items = "edgelist"
  )
edge_list <-
  edge_list %>%
  select(upia, upib, marker) %>%
  distinct()
bipart_graph <-
  edge_list %>%
  # Create tidy graph
  as_tbl_graph(directed = FALSE) %>%
  mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
attr(bipart_graph, "type") <- "bipartite"

test_that("CreateCellGraphObject works as expected", {
  cg <- CreateCellGraphObject(cellgraph = bipart_graph)
  expect_s4_class(cg, "CellGraph")
})

test_that("CellGraph initialize supplies default slot values", {
  cg <- methods::new("CellGraph")
  expect_s4_class(cg, "CellGraph")
  expect_null(cg@cellgraph)
  expect_null(cg@counts)
  expect_null(cg@layout)
  expect_identical(cg@layers, list())
  expect_equal(cg@meta.data, data.frame())
  expect_identical(cg@reductions, list())

  cg <- methods::new("CellGraph", cellgraph = bipart_graph)
  expect_identical(cg@layers, list())
  expect_identical(cg@reductions, list())
  expect_equal(rownames(cg@meta.data), bipart_graph %>% dplyr::pull(name))
})

test_that("AddMetaData fills per-node meta.data on empty tables", {
  nodes <- bipart_graph %>% dplyr::pull(name)
  cg <- CreateCellGraphObject(cellgraph = bipart_graph)
  slot(cg, "meta.data") <- data.frame()
  cg <- SeuratObject::AddMetaData(cg, metadata = seq_along(nodes), col.name = "idx")
  expect_equal(nrow(cg@meta.data), length(nodes))
  expect_equal(rownames(cg@meta.data), nodes)
  expect_equal(cg@meta.data$idx, seq_along(nodes))
})

test_that("CreateCellGraphObject accepts named layout lists", {
  layout <- tibble::tibble(x = seq_len(length(bipart_graph)))

  cg <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    layout = list(example_layout = layout)
  )

  expect_s3_class(cg@layout$example_layout, "data.frame")
  expect_false(inherits(cg@layout$example_layout, "tbl_df"))
  expect_equal(cg@layout$example_layout$x, layout$x)
  expect_equal(colnames(cg@layout$example_layout), "x")
  expect_equal(rownames(cg@layout$example_layout), bipart_graph %>% pull(name))
})

test_that("CreateCellGraphObject fails when invalid input is provided", {
  expect_error(CreateCellGraphObject(cellgraph = "Invalid input"))
  expect_error(CreateCellGraphObject(cellgraph = bipart_graph, counts = "Invalid input"))
  expect_error(CreateCellGraphObject(cellgraph = bipart_graph, layout = "Invalid input"))
  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      layout = list(example_layout = data.frame(x = 1))
    )
  )
  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      layout = list(tibble::tibble(x = 1))
    ),
    "must be named"
  )
  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      layers = list(counts = matrix(1))
    )
  )
})

node_names <- bipart_graph %>% dplyr::pull(name)
n_nodes <- length(node_names)

make_counts <- function(row_order = node_names) {
  mat <- Matrix::Matrix(
    seq_len(length(row_order) * 3),
    nrow = length(row_order),
    ncol = 3,
    sparse = TRUE
  )
  mat <- as(mat, "dgCMatrix")
  rownames(mat) <- row_order
  colnames(mat) <- c("m1", "m2", "m3")
  mat
}

test_that("CreateCellGraphObject aligns shuffled counts by node name", {
  shuffled_names <- rev(node_names)
  counts <- make_counts(shuffled_names)
  cg <- CreateCellGraphObject(cellgraph = bipart_graph, counts = counts)
  expect_equal(rownames(cg@counts), node_names)
  expect_equal(as.matrix(cg@counts), as.matrix(counts[node_names, ]))
})

test_that("CreateCellGraphObject aligns layouts with a name column", {
  layout <- tibble::tibble(
    name = rev(node_names),
    x = seq_len(n_nodes)
  )
  cg <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    layout = list(example_layout = layout)
  )
  expect_equal(cg@layout$example_layout$x, as.integer(match(node_names, rev(node_names))))
  expect_false("name" %in% colnames(cg@layout$example_layout))
  expect_equal(rownames(cg@layout$example_layout), node_names)
})

test_that("CreateCellGraphObject aligns layouts by row names", {
  layout <- data.frame(
    x = seq_len(n_nodes),
    y = seq_len(n_nodes),
    row.names = rev(node_names)
  )
  cg <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    layout = list(example_layout = layout)
  )
  expect_equal(rownames(cg@layout$example_layout), node_names)
  expect_equal(cg@layout$example_layout$x, layout[node_names, "x"])
})

test_that("CreateCellGraphObject stores layers, meta.data and reductions", {
  counts <- make_counts()
  layer_mat <- matrix(
    seq_len(n_nodes * 2),
    nrow = n_nodes,
    ncol = 2,
    dimnames = list(rev(node_names), c("f1", "f2"))
  )
  meta <- tibble::tibble(
    name = rev(node_names),
    cluster = rep(c("a", "b"), length.out = n_nodes),
    score = seq_len(n_nodes)
  )
  embeddings <- matrix(
    seq_len(n_nodes * 2),
    nrow = n_nodes,
    ncol = 2,
    dimnames = list(rev(node_names), NULL)
  )
  dr <- CreateNodeDimReducObject(
    embeddings = embeddings,
    key = "PC_",
    method = "pca",
    stdev = c(2, 1)
  )

  cg <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = counts,
    layers = list(data = layer_mat),
    meta.data = meta,
    reductions = list(pca = dr)
  )

  expect_equal(SeuratObject::Layers(cg), c("counts", "data"))
  expect_equal(rownames(SeuratObject::LayerData(cg, layer = "data")), node_names)
  expect_equal(rownames(cg@meta.data), node_names)
  expect_equal(cg@meta.data$cluster[1], meta$cluster[match(node_names[1], meta$name)])
  expect_equal(names(cg@reductions), "pca")
  expect_equal(names(CellGraphData(cg, slot = "meta_data")), names(cg@meta.data))
  expect_equal(rownames(SeuratObject::Embeddings(cg, reduction = "pca")), node_names)
  expect_s4_class(cg[["pca"]], "NodeDimReduc")
})

test_that("CellGraphData<- stores the provided value and remaps shuffled graphs", {
  cg <- CreateCellGraphObject(cellgraph = bipart_graph, counts = make_counts())

  # The replacement holds different values than the counts it replaces, so that
  # the setter has to store the value it was given
  replacement <- make_counts(rev(node_names)) * 10
  CellGraphData(cg, slot = "counts") <- replacement
  expect_equal(rownames(CellGraphData(cg, slot = "counts")), node_names)
  expect_equal(
    as.matrix(CellGraphData(cg, slot = "counts")),
    as.matrix(replacement[node_names, ])
  )

  shuffled_graph <- bipart_graph %N>% dplyr::arrange(dplyr::desc(name))
  attr(shuffled_graph, "type") <- "bipartite"
  CellGraphData(cg, slot = "cellgraph") <- shuffled_graph
  new_names <- shuffled_graph %>% dplyr::pull(name)
  expect_equal(rownames(cg@counts), new_names)
  expect_equal(as.matrix(cg@counts), as.matrix(replacement[new_names, ]))
})

test_that("LayerData and AddMetaData work on CellGraph objects", {
  cg <- CreateCellGraphObject(cellgraph = bipart_graph, counts = make_counts())
  new_layer <- matrix(
    1,
    nrow = n_nodes,
    ncol = 1,
    dimnames = list(node_names, "z")
  )
  SeuratObject::LayerData(cg, layer = "scaled") <- new_layer
  expect_equal(SeuratObject::Layers(cg), c("counts", "scaled"))
  expect_equal(as.vector(SeuratObject::LayerData(cg, layer = "scaled")), rep(1, n_nodes))

  cg <- SeuratObject::AddMetaData(cg, metadata = seq_len(n_nodes), col.name = "idx")
  expect_equal(cg@meta.data$idx, seq_len(n_nodes))

  cg[["umap"]] <- CreateNodeDimReducObject(
    embeddings = matrix(
      runif(n_nodes * 2),
      nrow = n_nodes,
      dimnames = list(node_names, NULL)
    ),
    key = "UMAP_",
    method = "umap"
  )
  expect_equal(names(cg@reductions), "umap")
})

test_that("matrix layers may share feature names", {
  counts <- make_counts()
  overlapping_layer <- matrix(
    100 + seq_len(n_nodes),
    ncol = 1,
    dimnames = list(node_names, "m1")
  )

  expect_no_error(
    cg <- CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      layers = list(data = overlapping_layer)
    )
  )
  expect_no_error(
    SeuratObject::LayerData(cg, layer = "scaled") <- overlapping_layer
  )
  expect_equal(colnames(cg@counts), c("m1", "m2", "m3"))
  expect_equal(colnames(cg@layers$data), "m1")
  expect_equal(colnames(cg@layers$scaled), "m1")
  expect_equal(
    SeuratObject::FetchData(cg, vars = "m1", layer = "counts")$m1,
    as.numeric(counts[, "m1"])
  )
  expect_equal(
    SeuratObject::FetchData(cg, vars = "m1", layer = "data")$m1,
    as.numeric(overlapping_layer[, "m1"])
  )
})

test_that("constructor rejects variable name collisions across data sources", {
  counts <- make_counts()
  duplicate_counts <- counts
  colnames(duplicate_counts)[2] <- colnames(duplicate_counts)[1]
  expect_error(
    CreateCellGraphObject(cellgraph = bipart_graph, counts = duplicate_counts),
    "Feature names in counts must be unique"
  )

  colliding_meta <- data.frame(
    m1 = seq_len(n_nodes),
    row.names = node_names
  )
  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      meta.data = colliding_meta
    ),
    "present in both meta.data and counts/layers"
  )

  colliding_graph <- bipart_graph %N>%
    dplyr::mutate(m1 = seq_len(n_nodes))
  expect_error(
    CreateCellGraphObject(cellgraph = colliding_graph, counts = counts),
    "present in both cellgraph node table and counts/layers"
  )

  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      meta.data = data.frame(
        node_type = rep("umi", n_nodes),
        row.names = node_names
      )
    ),
    "cellgraph node table and meta.data"
  )

  reduction <- CreateNodeDimReducObject(
    embeddings = matrix(
      seq_len(n_nodes),
      ncol = 1,
      dimnames = list(node_names, NULL)
    ),
    key = "m"
  )
  colnames(reduction@embeddings) <- "m1"
  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      reductions = list(pca = reduction)
    ),
    "reduction 'pca'.*counts/layers"
  )

  expect_error(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      reductions = list(first = reduction, second = reduction)
    ),
    "reduction 'first'.*reduction 'second'"
  )
})

test_that("CellGraph setters reject variable name collisions", {
  cg <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = make_counts(),
    meta.data = data.frame(cluster = rep("a", n_nodes), row.names = node_names)
  )

  expect_error(
    SeuratObject::AddMetaData(cg, seq_len(n_nodes), col.name = "m1"),
    "meta.data and counts/layers"
  )
  expect_no_error(
    cg <- SeuratObject::AddMetaData(cg, seq_len(n_nodes), col.name = "cluster")
  )

  colliding_layer <- matrix(
    seq_len(n_nodes),
    ncol = 1,
    dimnames = list(node_names, "cluster")
  )
  expect_error(
    SeuratObject::LayerData(cg, layer = "data") <- colliding_layer,
    "meta.data and counts/layers"
  )

  colliding_counts <- make_counts()
  colnames(colliding_counts)[1] <- "cluster"
  expect_error(
    CellGraphData(cg, slot = "counts") <- colliding_counts,
    "meta.data and counts/layers"
  )

  colliding_graph <- bipart_graph %N>%
    dplyr::mutate(m1 = seq_len(n_nodes))
  expect_error(
    CellGraphData(cg, slot = "cellgraph") <- colliding_graph,
    "cellgraph node table and counts/layers"
  )

  reduction <- CreateNodeDimReducObject(
    embeddings = matrix(
      seq_len(n_nodes),
      ncol = 1,
      dimnames = list(node_names, NULL)
    ),
    key = "m"
  )
  colnames(reduction@embeddings) <- "m1"
  expect_error(
    cg[["pca"]] <- reduction,
    "reduction 'pca'.*counts/layers"
  )

  expect_error(
    CellGraphData(cg, slot = "layers") <- list(data = colliding_layer),
    "meta.data and counts/layers"
  )
  expect_error(
    CellGraphData(cg, slot = "meta.data") <- data.frame(
      m1 = seq_len(n_nodes),
      row.names = node_names
    ),
    "meta.data and counts/layers"
  )
  expect_error(
    CellGraphData(cg, slot = "reductions") <- list(pca = reduction),
    "reduction 'pca'.*counts/layers"
  )
})

test_that("subset.CellGraph keeps node-level slots aligned", {
  counts <- make_counts()
  layout <- tibble::tibble(x = seq_len(n_nodes), y = seq_len(n_nodes))
  layer_mat <- matrix(
    seq_len(n_nodes),
    ncol = 1,
    dimnames = list(node_names, "f1")
  )
  cg <- CreateCellGraphObject(
    cellgraph = bipart_graph,
    counts = counts,
    layout = list(xy = layout),
    layers = list(data = layer_mat),
    meta.data = data.frame(grp = rep("a", n_nodes), row.names = node_names)
  )
  keep <- node_names[seq_len(min(50, n_nodes))]
  cg_small <- subset(cg, nodes = keep)
  small_names <- cg_small@cellgraph %>% dplyr::pull(name)
  expect_equal(rownames(cg_small@counts), small_names)
  expect_equal(rownames(cg_small@layout$xy), small_names)
  expect_equal(rownames(cg_small@layers$data), small_names)
  expect_equal(rownames(cg_small@meta.data), small_names)
  expect_true(all(small_names %in% keep))
  expect_equal(cg_small@layout$xy$x, layout$x[match(small_names, node_names)])
})

test_that("subset.CellGraph accepts a single node", {
  cg <- CreateCellGraphObject(cellgraph = bipart_graph)
  one <- node_names[1]
  cg_one <- subset(cg, nodes = one)
  expect_equal(igraph::gorder(cg_one@cellgraph), 1)
  expect_equal(cg_one@cellgraph %>% dplyr::pull(name), one)
})

test_that("KeepLargestComponent.CellGraph accepts a single-node component", {
  g <- tidygraph::tbl_graph(
    nodes = data.frame(name = c("a", "b", "c"), stringsAsFactors = FALSE),
    edges = data.frame(from = integer(), to = integer())
  )
  attr(g, "type") <- "single"
  cg <- CreateCellGraphObject(cellgraph = g)
  expect_no_error(cg_largest <- KeepLargestComponent(cg, verbose = FALSE))
  expect_equal(igraph::gorder(cg_largest@cellgraph), 1)
})

test_that("legacy layout tables without row names get node IDs on access", {
  layout <- tibble::tibble(x = seq_len(n_nodes), y = seq_len(n_nodes))
  cg <- CreateCellGraphObject(cellgraph = bipart_graph, layout = list(xy = layout))
  cg@layout$xy <- layout

  cg_small <- subset(cg, nodes = node_names[1:3])
  expect_equal(rownames(cg_small@layout$xy), node_names[1:3])
})
