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

test_that("CreateCellGraphObject fails when invalid input is provided", {
  expect_error(CreateCellGraphObject(cellgraph = "Invalid input"))
  expect_error(CreateCellGraphObject(cellgraph = bipart_graph, counts = "Invalid input"))
  expect_error(CreateCellGraphObject(cellgraph = bipart_graph, layout = "Invalid input"))
})
