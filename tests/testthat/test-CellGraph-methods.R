library(tidygraph)

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

cg <- CreateCellGraphObject(cellgraph = bipart_graph)

# CellGraphData method
test_that("CellGraphData works as expected", {
  expect_no_error(cellgraph <- CellGraphData(cg))
  expect_s3_class(cellgraph, class = "tbl_graph")
  expect_equal(cellgraph %>% igraph::gsize(), 68255)
  expect_equal(cellgraph %>% length(), 16800)
})

test_that("CellGraphData fails when invalid input is provided", {
  expect_error(CellGraphData("Invalid input"))
  expect_error(CellGraphData(cg, slot = "Invalid"))
})

# CellGraphData<- method
test_that("CellGraphData<- works as expected", {
  expect_no_error(CellGraphData(cg) <- CellGraphData(cg))
  expect_no_error(cellgraph <- CellGraphData(cg))
  expect_s3_class(cellgraph, class = "tbl_graph")
  expect_equal(cellgraph %>% igraph::gsize(), 68255)
  expect_equal(cellgraph %>% length(), 16800)
})

test_that("CellGraphData<- fails when invalid input is provided", {
  expect_error(CellGraphData(cg) <- "Invalid")
  expect_error(CellGraphData(cg, slot = "counts") <- "Invalid")
  expect_error(CellGraphData(cg, slot = "layout") <- "Invalid")
})

# show method
test_that("show.CellGraph works as expected", {
  msg <- capture_output(show(cg))
  expect_equal(msg, "A CellGraph object containing a bipartite graph with 16800 nodes and 68255 edges")
})
