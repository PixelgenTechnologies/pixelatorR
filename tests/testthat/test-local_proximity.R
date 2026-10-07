library(dplyr)
se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
  LoadCellGraphs(cells = colnames(.)[1])
cg <- CellGraphs(se)[[1]]
set.seed(123)
nodes <- cg@cellgraph %>%
  pull(name) %>%
  sample(20000)
cg_small <- subset(cg, nodes = nodes)

test_that("local_proximity works as expected", {
  # Single marker
  expect_no_error({
    score <- local_proximity(
      object = cg_small,
      markers = "B2M"
    )
  })

  expect_equal(
    score %>% sort(decreasing = TRUE) %>% head(2),
    c(
      `21556371155279213` = 2.05254858649264,
      `62286657368003726` = 2.03216767866132
    )
  )

  # Two markers
  expect_no_error({
    score <- local_proximity(
      object = cg_small,
      markers = c("B2M", "HLA-ABC")
    )
  })

  expect_equal(
    score %>% sort(decreasing = TRUE) %>% head(2),
    c(
      `8863920485090395` = 1.9603110060193,
      `5291964611270309` = 1.82087790932478
    )
  )

  # Several markers
  expect_no_error({
    score <- local_proximity(
      object = cg_small,
      markers = c("B2M", "HLA-ABC", "CD45")
    )
  })

  expect_equal(
    score %>% sort(decreasing = TRUE) %>% head(2),
    c(
      `68217662764440214` = 1.16803078459537,
      `68426598618843577` = 1.1496758370931
    )
  )

  # any mode
  expect_no_error({
    score <- local_proximity(
      object = cg_small,
      markers = c("B2M", "HLA-ABC"),
      mode = "any"
    )
  })

  expect_equal(
    score %>% sort(decreasing = TRUE) %>% head(2),
    c(
      `62286657368003726` = 2.45069970675555,
      `8863920485090395` = 2.3533575067355
    )
  )

  # self-clustering mode
  expect_no_error({
    score <- local_proximity(
      object = cg_small,
      markers = c("B2M", "HLA-ABC"),
      mode = "self-clustering"
    )
  })

  expect_type(score, "double")
  expect_true(all(dim(score) == c(20000, 2)))

  # permutations
  expect_no_error({
    score <- local_proximity(
      object = cg_small,
      markers = c("B2M", "HLA-ABC"),
      method = "permutation",
      iterations = 10
    )
  })

  expect_equal(
    score %>% sort(decreasing = TRUE) %>% head(2),
    c(
      `4819662606092360` = 2.15200309344505,
      `62286657368003726` = 2.09953567355091
    )
  )
})

test_that("local_proximity fails with invalid input", {
  expect_error(local_proximity(cg_small, markers = "Invalid"))
  expect_error(local_proximity(cg_small, markers = "B2M", method = "Invalid"))
  expect_error(local_proximity(cg_small, markers = "B2M", mode = "Invalid"))
  expect_error(local_proximity(cg_small, markers = "B2M", iterations = "Invalid"))
  expect_error(local_proximity(cg_small, markers = "B2M", k = "Invalid"))
  expect_error(local_proximity(cg_small, markers = "B2M", seed = "Invalid"))
  expect_error(local_proximity(cg_small, markers = "B2M", method = "permutation", mode = "self-clustering"))
})
