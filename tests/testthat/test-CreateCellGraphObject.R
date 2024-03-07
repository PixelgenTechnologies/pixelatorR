edge_list <-
  ReadMPX_item(
    system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR"),
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

test_that("CreateCellGraphObject fails when invalid input is provided", {
  expect_error(CreateCellGraphObject(cellgraph = "Invalid input"), "'cellgraph' must be a non-empty")
  expect_error(CreateCellGraphObject(cellgraph = bipart_graph, counts = "Invalid input"), "'counts' must be a non-empty 'dgCMatrix' object")
  expect_error(CreateCellGraphObject(cellgraph = bipart_graph, layout = "Invalid input"), "'layout' must be a non-empty 'tbl_df' object")
})
