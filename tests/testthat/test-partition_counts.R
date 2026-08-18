library(dplyr)
library(tidygraph)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE)
cg <- CellGraphs(se)[[4]]
partition <- rep(c("A", "B"), length.out = length(cg@cellgraph))
cg@cellgraph <- cg@cellgraph %N>%
  mutate(partition = partition)

expected_counts <- structure(
  c(8516, 8588, 1626, 1644, 265, 307),
  dim = 2:3,
  dimnames = list(c("A", "B"), c("CD44", "CD3e", "CD4"))
)

test_that("partition_counts works as expected", {
  expect_no_error(
    res <- partition_counts(cg, partition = partition)
  )

  expect_equal(
    as.matrix(res[, colnames(expected_counts)]),
    expected_counts
  )

  expect_no_error(
    res <- partition_counts(cg, partition_column = "partition")
  )

  expect_equal(
    as.matrix(res[, colnames(expected_counts)]),
    expected_counts
  )
})


test_that("partition_counts fails with invalid input", {
  expect_error(
    partition_counts("Invalid")
  )

  expect_error(
    partition_counts(cg, partition = 1:10)
  )

  expect_error(
    partition_counts(cg, partition = c("A", "B"))
  )

  expect_error(
    partition_counts(cg, partition_column = "Invalid")
  )
})
