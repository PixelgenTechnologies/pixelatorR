options(Seurat.object.assay.version = "v5")
library(dplyr)
library(tidygraph)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[4], verbose = FALSE)
cg <- CellGraphs(se)[[4]]

test_that("partition_counts works as expected with partition vector", {
  nodes <- cg@cellgraph %N>% pull(name)
  partition <- rep(c("A", "B"), length.out = length(nodes))

  expect_no_error(counts <- partition_counts(cg, partition = partition))

  expect_equal(rownames(counts), c("A", "B"))
  expect_equal(colnames(counts), colnames(cg@counts))
  expect_equal(
    as.matrix(counts[, c("CD44", "CD3e", "CD4")]),
    structure(
      c(8516, 8588, 1626, 1644, 265, 307),
      dim = 2:3,
      dimnames = list(c("A", "B"), c("CD44", "CD3e", "CD4"))
    )
  )
  expect_equal(
    as.numeric(Matrix::colSums(counts)),
    as.numeric(Matrix::colSums(cg@counts))
  )
})

test_that("partition_counts works as expected with partition_column", {
  cg_labeled <- cg
  cg_labeled@cellgraph <- cg_labeled@cellgraph %N>%
    mutate(group = rep(c("group1", "group2"), length.out = length(cg_labeled@cellgraph)))

  expect_no_error(
    counts <- partition_counts(cg_labeled, partition_column = "group")
  )
  expect_equal(rownames(counts), c("group1", "group2"))
  expect_equal(
    as.numeric(Matrix::colSums(counts)),
    as.numeric(Matrix::colSums(cg_labeled@counts))
  )
})

test_that("partition_counts fails with invalid input", {
  nodes <- cg@cellgraph %N>% pull(name)

  expect_error(partition_counts("Invalid"))
  expect_error(partition_counts(cg))
  expect_error(
    partition_counts(cg, partition = rep("A", length(nodes)), partition_column = "group")
  )
  expect_error(partition_counts(cg, partition = "Invalid"))
  expect_error(partition_counts(cg, partition_column = "Invalid"))
})
