library(tidygraph)
set.seed(123)
cg <- simulate_bipartite_graph(n_nodes = 1e3)
cg <- add_binary_marker_counts(cg, markers = paste0("marker", 1:2))
cg@cellgraph <- cg@cellgraph %N>%
  mutate(partition = sample(c("A", "B"), n(), replace = TRUE))
partition <- cg@cellgraph %N>% pull(partition)

test_that("partition_counts works as expected", {
  expect_no_error(
    res <- partition_counts(cg, partition = partition)
  )

  expect_equal(
    new(
      "dgCMatrix",
      i = c(0L, 1L, 0L, 1L),
      p = c(0L, 2L, 4L),
      Dim = c(2L, 2L),
      Dimnames = list(c("A", "B"), c("marker1", "marker2")),
      x = c(265, 275, 257, 227),
      factors = list()
    ),
    res
  )

  expect_no_error(
    res <- partition_counts(cg, partition_column = "partition")
  )

  expect_equal(
    new(
      "dgCMatrix",
      i = c(0L, 1L, 0L, 1L),
      p = c(0L, 2L, 4L),
      Dim = c(2L, 2L),
      Dimnames = list(c("A", "B"), c("marker1", "marker2")),
      x = c(265, 275, 257, 227),
      factors = list()
    ),
    res
  )
})


test_that("partition_counts fails with invalid input", {
  expect_error(
    res <- partition_counts("Invalid")
  )

  expect_error(
    res <- partition_counts(cg, partition = 1:10)
  )

  expect_error(
    res <- partition_counts(cg, partition = c("A", "B"))
  )

  expect_error(
    res <- partition_counts(cg, partition_column = "Invalid")
  )
})
