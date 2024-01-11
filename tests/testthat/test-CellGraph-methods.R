edge_list <-
  ReadMPX_item(
    system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR"),
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

cg <- CreateCellGraphObject(cellgraph = bipart_graph)

# GetCellGraphData method
test_that("GetCellGraphData works as expected", {
  cellgraph <- GetCellGraphData(cg)
  expect_s3_class(cellgraph, class = "tbl_graph")
  expect_equal(cellgraph %>% gsize(), 46655)
  expect_equal(cellgraph %>% length(), 26964)
})

test_that("GetCellGraphData fails when invalid input is provided", {
  expect_error(GetCellGraphData("Invalid input"), "Invalid class character")
  expect_error(GetCellGraphData(cg, slot = "Invalid"), "slot must be one of cellgraph, counts, layout")
})

# show method
test_that("show.CellGraph works as expected", {
  msg <- capture_output(show(cg))
  expect_equal(msg, "A CellGraph object containing a bipartite graph with 26964 nodes and 46655 edges ")
})
