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
      `11003117224941247-umi2` = 2.92922642517793,
      `67262613737613740-umi1` = 2.66192180611334
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
      `67262613737613740-umi1` = 2.8073549220576,
      `8936220015879929-umi2` = 2.4594316186373
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
      `34669962850295590-umi1` = 0.584962500721156,
      `61322123765323003-umi1` = 0.584962500721156
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
      `57061743620873557-umi1` = 3.08746284125034,
      `40552551807181406-umi1` = 3.08746284125034
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
      `67262613737613740-umi1` = 2.8073549220576,
      `8936220015879929-umi2` = 2.4594316186373
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
