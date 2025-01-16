pxl_file <- system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")
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

test_that("CreateCellGraphAssay5 works as expected", {
  expect_no_error({
    cg_assay <- CreateCellGraphAssay5(counts = mat, cellgraphs = bipartite_graphs)
  })
  expect_s4_class(cg_assay, "CellGraphAssay5")
  expect_no_error({
    cg_assay <- CreateCellGraphAssay5(counts = mat, cellgraphs = bipartite_graphs)
  })
  expect_s4_class(cg_assay, "CellGraphAssay5")
})

test_that("CreateCellGraphAssay5 fails when invalid input is provided", {
  expect_error(CreateCellGraphAssay5(counts = "Invalid input", cellgraphs = bipartite_graphs))
  expect_error(CreateCellGraphAssay5(counts = mat, cellgraphs = "Invalid input"))
})
