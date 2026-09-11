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

# cli wraps abort messages at the console width, so collapse whitespace before
# matching to keep the expected text independent of where the line breaks fall.
expect_error_text <- function(object, pattern) {
  message <- tryCatch(
    {
      object
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (is.null(message)) {
    fail("Expected an error, but none was thrown.")
  }
  expect_match(gsub("[[:space:]]+", " ", message), pattern)
}

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
  expect_equal(cg@nodes, character())

  cg <- methods::new("CellGraph", cellgraph = bipart_graph)
  expect_identical(cg@layers, list())
  expect_identical(cg@reductions, list())
  expect_equal(cg@nodes, bipart_graph %>% dplyr::pull(name))
  expect_equal(nrow(cg@meta.data), length(cg@nodes))
  expect_lt(.row_names_info(cg@meta.data), 0L)
})

test_that("AddMetaData fills per-node meta.data on empty tables", {
  nodes <- bipart_graph %>% dplyr::pull(name)
  cg <- CreateCellGraphObject(cellgraph = bipart_graph)
  slot(cg, "meta.data") <- data.frame()
  cg <- SeuratObject::AddMetaData(cg, metadata = seq_along(nodes), col.name = "idx")
  expect_equal(nrow(cg@meta.data), length(nodes))
  expect_equal(SeuratObject::Cells(cg), nodes)
  expect_lt(.row_names_info(cg@meta.data), 0L)
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
  expect_equal(SeuratObject::Cells(cg), bipart_graph %>% pull(name))
  expect_lt(.row_names_info(cg@layout$example_layout), 0L)
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
  expect_null(rownames(cg@counts))
  expect_equal(SeuratObject::Cells(cg), node_names)
  expect_equal(as.matrix(cg@counts), {
    expected <- as.matrix(counts[node_names, ])
    dimnames(expected) <- dimnames(as.matrix(cg@counts))
    expected
  })
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
  expect_lt(.row_names_info(cg@layout$example_layout), 0L)
  expect_equal(SeuratObject::Cells(cg), node_names)
})

test_that("CreateCellGraphObject aligns MPX layouts without A/B suffixes", {
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

  layout_named <- tibble::tibble(
    name = c("umi2", "umi1"),
    x = c(20, 10),
    y = c(2, 1)
  )
  cg <- CreateCellGraphObject(
    cellgraph = g,
    layout = list(pmds_3d = layout_named)
  )
  expect_equal(cg@layout$pmds_3d$x, c(10, 10, 20, 20))
  expect_equal(cg@layout$pmds_3d$y, c(1, 1, 2, 2))
  expect_false("name" %in% colnames(cg@layout$pmds_3d))
  expect_equal(SeuratObject::Cells(cg), suffixed)

  layout_rownames <- data.frame(
    x = c(10, 20),
    y = c(1, 2),
    row.names = c("umi1", "umi2")
  )
  CellGraphData(cg, slot = "layout") <- list(pmds_3d = layout_rownames)
  expect_equal(cg@layout$pmds_3d$x, c(10, 10, 20, 20))
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
  expect_lt(.row_names_info(cg@layout$example_layout), 0L)
  expect_equal(SeuratObject::Cells(cg), node_names)
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
  expect_null(rownames(cg@layers$data))
  expect_lt(.row_names_info(cg@meta.data), 0L)
  expect_equal(SeuratObject::Cells(cg), node_names)
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
  expect_null(rownames(CellGraphData(cg, slot = "counts")))
  expect_equal(SeuratObject::Cells(cg), node_names)
  expect_equal(
    as.numeric(as.matrix(CellGraphData(cg, slot = "counts"))),
    as.numeric(as.matrix(replacement[node_names, ]))
  )

  shuffled_graph <- bipart_graph %N>% dplyr::arrange(dplyr::desc(name))
  attr(shuffled_graph, "type") <- "bipartite"
  CellGraphData(cg, slot = "cellgraph") <- shuffled_graph
  new_names <- shuffled_graph %>% dplyr::pull(name)
  expect_null(rownames(cg@counts))
  expect_equal(SeuratObject::Cells(cg), new_names)
  expect_equal(as.numeric(as.matrix(cg@counts)), as.numeric(as.matrix(replacement[new_names, ])))
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

test_that("AddMetaData annotates a subset of nodes", {
  cg <- CreateCellGraphObject(cellgraph = bipart_graph)

  partial <- c(1, 3)
  names(partial) <- node_names[c(1, 3)]
  cg <- SeuratObject::AddMetaData(cg, metadata = partial, col.name = "score")
  expect_equal(nrow(cg@meta.data), n_nodes)
  expect_equal(SeuratObject::Cells(cg), node_names)
  expect_equal(cg@meta.data$score[c(1, 3)], c(1, 3))
  expect_true(all(is.na(cg@meta.data$score[-c(1, 3)])))

  # Names that are not graph nodes are dropped
  partial_df <- data.frame(
    group = c("a", "b"),
    row.names = c(node_names[2], "not_a_node")
  )
  cg <- SeuratObject::AddMetaData(cg, metadata = partial_df)
  expect_equal(cg@meta.data$group[2], "a")
  expect_true(all(is.na(cg@meta.data$group[-2])))

  expect_error(
    SeuratObject::AddMetaData(
      cg,
      metadata = data.frame(x = 1, row.names = "not_a_node")
    ),
    "No node"
  )
  expect_error(
    SeuratObject::AddMetaData(
      cg,
      metadata = c(a = 1, a = 2),
      col.name = "dup"
    ),
    "must be unique"
  )
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
  expect_error_text(
    CreateCellGraphObject(cellgraph = bipart_graph, counts = duplicate_counts),
    "Feature names in counts must be unique"
  )

  colliding_meta <- data.frame(
    m1 = seq_len(n_nodes),
    row.names = node_names
  )
  expect_error_text(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      meta.data = colliding_meta
    ),
    "present in both meta.data and counts/layers"
  )

  colliding_graph <- bipart_graph %N>%
    dplyr::mutate(m1 = seq_len(n_nodes))
  expect_error_text(
    CreateCellGraphObject(cellgraph = colliding_graph, counts = counts),
    "present in both cellgraph node table and counts/layers"
  )

  expect_error_text(
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
  expect_error_text(
    CreateCellGraphObject(
      cellgraph = bipart_graph,
      counts = counts,
      reductions = list(pca = reduction)
    ),
    "reduction 'pca'.*counts/layers"
  )

  expect_error_text(
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

  expect_error_text(
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
  expect_error_text(
    SeuratObject::LayerData(cg, layer = "data") <- colliding_layer,
    "meta.data and counts/layers"
  )

  colliding_counts <- make_counts()
  colnames(colliding_counts)[1] <- "cluster"
  expect_error_text(
    CellGraphData(cg, slot = "counts") <- colliding_counts,
    "meta.data and counts/layers"
  )

  colliding_graph <- bipart_graph %N>%
    dplyr::mutate(m1 = seq_len(n_nodes))
  expect_error_text(
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
  expect_error_text(
    cg[["pca"]] <- reduction,
    "reduction 'pca'.*counts/layers"
  )

  expect_error_text(
    CellGraphData(cg, slot = "layers") <- list(data = colliding_layer),
    "meta.data and counts/layers"
  )
  expect_error_text(
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
  expect_equal(SeuratObject::Cells(cg_small), small_names)
  expect_null(rownames(cg_small@counts))
  expect_null(rownames(cg_small@layers$data))
  expect_lt(.row_names_info(cg_small@layout$xy), 0L)
  expect_lt(.row_names_info(cg_small@meta.data), 0L)
  expect_equal(nrow(cg_small@counts), length(small_names))
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

test_that("layout tables without row names get node IDs on subset", {
  layout <- tibble::tibble(x = seq_len(n_nodes), y = seq_len(n_nodes))
  cg <- CreateCellGraphObject(cellgraph = bipart_graph, layout = list(xy = layout))
  cg@layout$xy <- layout

  cg_small <- subset(cg, nodes = node_names[1:3])
  expect_equal(SeuratObject::Cells(cg_small), node_names[1:3])
  expect_equal(cg_small@layout$xy$x, layout$x[1:3])
})

make_single_graph <- function(node_ids = NULL) {
  n <- if (is.null(node_ids)) 3L else length(node_ids)
  nodes <- data.frame(node_type = rep("A", n), stringsAsFactors = FALSE)
  if (!is.null(node_ids)) {
    nodes$name <- node_ids
  }
  g <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = data.frame(from = seq_len(n - 1L), to = seq_len(n - 1L) + 1L)
  )
  attr(g, "type") <- "single"
  g
}

test_that("subset.CellGraph keeps sequential IDs on graphs without a name attribute", {
  cg <- CreateCellGraphObject(cellgraph = make_single_graph())
  expect_equal(cg@cellgraph %>% dplyr::pull(name), c("1", "2", "3"))

  cg_small <- subset(cg, nodes = c("2", "3"))
  expect_equal(cg_small@cellgraph %>% dplyr::pull(name), c("2", "3"))
  expect_equal(SeuratObject::Cells(cg_small), c("2", "3"))
  expect_equal(nrow(cg_small@meta.data), 2)
})

test_that("KeepLargestComponent.CellGraph works for graphs without a name attribute", {
  g <- tidygraph::tbl_graph(
    nodes = data.frame(node_type = c("A", "A", "A"), stringsAsFactors = FALSE),
    edges = data.frame(from = 1L, to = 2L)
  )
  attr(g, "type") <- "single"
  cg <- CreateCellGraphObject(cellgraph = g)
  cg_largest <- KeepLargestComponent(cg, verbose = FALSE)
  expect_equal(sort(cg_largest@cellgraph %>% dplyr::pull(name)), c("1", "2"))
})

test_that("subset.CellGraph rebuilds empty meta.data with automatic rownames", {
  cg <- CreateCellGraphObject(cellgraph = make_single_graph(c("1", "2", "3")))
  auto_meta <- data.frame(x = 1:3)
  auto_meta$x <- NULL
  expect_true(.row_names_info(auto_meta) < 0L)
  cg@meta.data <- auto_meta

  cg_small <- subset(cg, nodes = c("2", "3"))
  expect_equal(SeuratObject::Cells(cg_small), c("2", "3"))
  expect_equal(nrow(cg_small@meta.data), 2)
})

test_that("subset.CellGraph keeps character integer node names", {
  cg <- CreateCellGraphObject(cellgraph = make_single_graph(c("1", "2", "3")))
  cg_small <- subset(cg, nodes = c("1", "3"))
  expect_equal(cg_small@cellgraph %>% dplyr::pull(name), c("1", "3"))
  expect_equal(SeuratObject::Cells(cg_small), c("1", "3"))
})

test_that("automatic layout rownames follow graph order, not 1:n IDs", {
  names <- c("umiA", "umiB", "umiC")
  layout <- data.frame(x = c(10, 20, 30), y = c(1, 2, 3))
  expect_true(.row_names_info(layout) < 0L)

  cg <- CreateCellGraphObject(
    cellgraph = make_single_graph(names),
    layout = list(xy = layout)
  )
  expect_lt(.row_names_info(cg@layout$xy), 0L)
  expect_equal(SeuratObject::Cells(cg), names)
  expect_equal(cg@layout$xy$x, c(10, 20, 30))
})

test_that("CellGraph objects from older versions report why they fail", {
  cg <- CreateCellGraphObject(cellgraph = bipart_graph, counts = make_counts())

  # Objects saved before the class gained layers, meta.data, and reductions
  # keep only the three original slots when they are read back from an RDS
  legacy <- cg
  attr(legacy, "layers") <- NULL
  attr(legacy, "meta.data") <- NULL
  attr(legacy, "reductions") <- NULL

  expect_error_text(print(legacy), "no layers, meta.data, and reductions slots")
  expect_error_text(print(legacy), "saved by pixelatorR 0.21.0 or earlier")
  expect_error_text(print(legacy), "pixelatorR@v0.20.1")

  expect_error_text(SeuratObject::Layers(legacy), "pixelatorR 0.21.0 or earlier")
  expect_error_text(SeuratObject::Cells(legacy), "pixelatorR 0.21.0 or earlier")
  expect_error_text(CellGraphData(legacy, slot = "counts"), "pixelatorR 0.21.0 or earlier")
  expect_error_text(SeuratObject::FetchData(legacy, vars = "m1"), "pixelatorR 0.21.0 or earlier")
  expect_error_text(subset(legacy, nodes = node_names[1]), "pixelatorR 0.21.0 or earlier")
  expect_error_text(
    SeuratObject::AddMetaData(legacy, metadata = seq_len(n_nodes), col.name = "idx"),
    "pixelatorR 0.21.0 or earlier"
  )
  expect_error_text(FetchLayoutData(legacy), "pixelatorR 0.21.0 or earlier")

  # Only the missing slots are named
  partial <- cg
  attr(partial, "reductions") <- NULL
  expect_error_text(print(partial), "no reductions slot\\b")
})
