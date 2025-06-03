library(pixelatorR)
se <- ReadMPX_Seurat(minimal_mpx_pxl_file(), verbose = FALSE)
se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE) %>%
  ComputeLayout(layout_method = "pmds", dim = 3, pivots = 50)
cg <- CellGraphs(se)[[1]]

# CellGraphData method
test_that("CellGraphData works as expected", {
  expect_no_error(cellgraph <- CellGraphData(cg))
  expect_s3_class(cellgraph, class = "tbl_graph")
  expect_equal(cellgraph %>% igraph::gsize(), 5138)
  expect_equal(cellgraph %>% length(), 2470)
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
  expect_equal(cellgraph %>% igraph::gsize(), 5138)
  expect_equal(cellgraph %>% length(), 2470)
})

test_that("CellGraphData<- fails when invalid input is provided", {
  expect_error(CellGraphData(cg) <- "Invalid")
  expect_error(CellGraphData(cg, slot = "counts") <- "Invalid")
  expect_error(CellGraphData(cg, slot = "layout") <- "Invalid")
})

# show method
test_that("show.CellGraph works as expected", {
  msg <- capture_output(show(cg))
  expect_equal(msg, "A CellGraph object containing a bipartite graph with 2470 nodes and 5138 edges\nNumber of markers:  79 \nLayouts: pmds_3d ")
})

# subset method
test_that("subset.CellGraph works as expected", {
  expect_no_error(cg_small <- subset(cg, nodes = rownames(cg@counts)[1:2000]))
  expect_equal(cg_small@layout$pmds_3d %>% dim(), c(2000, 3))
  expect_equal(cg_small@counts %>% dim(), c(2000, 79))
  expect_equal(cg_small@cellgraph %>% length(), 2000)
  expect_equal(cg_small@cellgraph %>% igraph::gsize(), 4227)
})
