edge_list <-
  ReadMPX_edgelist(
    system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR")
  )

test_that("edgelist_to_simple_Anode_graph works as expected", {
  edge_list_a_node <- edge_list %>%
    as_tibble() %>%
    edgelist_to_simple_Anode_graph()
  expect_type(edge_list_a_node, "list")
  expect_equal(edge_list_a_node %>% length(), 10)
  expect_equal(sapply(edge_list_a_node, function(x) class(x)[1]) %>% unique(), "tbl_graph")
})

test_that("edgelist_to_simple_Anode_graph fails when invalid input is provided", {
  expect_error(edgelist_to_simple_Anode_graph(object = tibble()), "edgelist must be a non-empty object")
  expect_error(edgelist_to_simple_Anode_graph(object = data.frame(x = seq(1, 10))), "One or several of 'upia', 'upib' are missing from edgelist")
})
