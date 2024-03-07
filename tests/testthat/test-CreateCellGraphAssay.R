pxl_file <- system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR")
mat <- ReadMPX_counts(pxl_file)
edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
components <- colnames(mat)
edgelist <-
  edgelist %>%
  select(upia, upib, marker, component) %>%
  distinct() %>%
  group_by(component) %>%
  group_split() %>%
  setNames(nm = components)
bipartite_graphs <- lapply(edgelist, function(x) {
  g <- x %>% tidygraph::as_tbl_graph(directed = FALSE)
  g <- g %>% mutate(node_type = case_when(name %in% x$upia ~ "A", TRUE ~ "B"))
  attr(g, "type") <- "bipartite"
  CreateCellGraphObject(cellgraph = g)
})

test_that("CreateCellGraphAssay works as expected", {
  expect_no_error({cg_assay <- CreateCellGraphAssay(counts = mat, cellgraphs = bipartite_graphs, overwrite = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_no_error({cg_assay <- CreateCellGraphAssay(counts = mat, cellgraphs = bipartite_graphs, arrow_dir = pxl_file, overwrite = TRUE)})
  expect_s4_class(cg_assay, "CellGraphAssay")
  expect_true(!is.na(cg_assay@arrow_dir))
  expect_equal(class(cg_assay@arrow_data), c("FileSystemDataset", "Dataset", "ArrowObject", "R6"))
})

test_that("CreateCellGraphAssay fails when invalid input is provided", {
  expect_error(CreateCellGraphAssay(counts = "Invalid input", cellgraphs = bipartite_graphs), "'counts' must be a matrix-like object")
  expect_error(CreateCellGraphAssay(counts = mat, cellgraphs = "Invalid input"), "'cellgraphs' must be a 'list'")
})

