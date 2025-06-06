edge_list <-
  ReadMPX_item(
    system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
    items = "edgelist"
  ) %>%
  filter(component == "RCVCMP0000217")

test_that("node_markers_counts works as expected", {
  # Test default settings
  node_marker_counts <-
    edge_list %>%
    node_markers_counts()
  expect_type(node_marker_counts, "integer")
  expect_equal(ncol(node_marker_counts), 79)
  expect_equal(nrow(node_marker_counts), 1395)

  # Test with k > 0
  node_marker_counts_k10 <-
    edge_list %>%
    node_markers_counts(k = 2)
  expect_type(node_marker_counts, "integer")
  expect_equal(ncol(node_marker_counts), 79)
  expect_equal(nrow(node_marker_counts), 1395)
})

test_that("node_markers_counts fails when invalid input is provided", {
  expect_error(node_markers_counts(component_edge_list = tibble()))
  expect_error(node_markers_counts(component_edge_list = data.frame(x = seq(1, 10))))
  expect_error(node_markers_counts(component_edge_list = tibble(x = seq(1, 10))))
})
