options(Seurat.object.assay.version = "v5")
library(dplyr)
library(tidygraph)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[4], verbose = FALSE)
cg <- CellGraphs(se)[[4]]

test_that("distance_from_node_set works as expected", {
  start_set <- cg@cellgraph %N>%
    pull(name) %>%
    head(1)

  expect_no_error(out <- distance_from_node_set(cg, start_set))

  node_data <- out@cellgraph %N>% as_tibble()
  expect_true("distance_from_seed" %in% names(node_data))
  expect_equal(node_data$distance_from_seed[match(start_set, node_data$name)], 0L)
  expect_equal(
    node_data$distance_from_seed[1:10],
    c(0L, 1L, 10L, 11L, 12L, 11L, 8L, 9L, 10L, 11L)
  )
  expect_equal(max(node_data$distance_from_seed, na.rm = TRUE), 21L)
  expect_equal(sum(is.na(node_data$distance_from_seed)), 0L)

  # Re-running should replace the existing distance column without error
  expect_no_error(distance_from_node_set(out, start_set, max_iter = 5L))
})

test_that("distance_from_node_set fails with invalid input", {
  start_set <- cg@cellgraph %N>%
    pull(name) %>%
    head(1)

  expect_error(distance_from_node_set("Invalid", start_set))
  expect_error(distance_from_node_set(cg, seed_nodes = "Invalid"))
  expect_error(distance_from_node_set(cg, seed_nodes = start_set, max_iter = "Invalid"))
  expect_error(distance_from_node_set(cg, seed_nodes = start_set, verbose = "Invalid"))
  expect_error(distance_from_node_set(cg, seed_nodes = character(0)))
})
