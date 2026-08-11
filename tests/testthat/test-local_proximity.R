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
      `11003117224941247-umi2` = 2.13410426952319,
      `9816302186898773-umi2` = 2.03158358758846
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
      `67262613737613740-umi1` = 2.10825438436169,
      `8936220015879929-umi2` = 1.94611096348858
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
      `41752602342049118-umi2` = 1.1737434287477,
      `61322123765323003-umi1` = 1.13209860331366
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
      `67262613737613740-umi1` = 2.43368047619204,
      `57061743620873557-umi1` = 2.27950747300082
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
      `67262613737613740-umi1` = 2.41503749927884,
      `11915219229172862-umi2` = 2.13750352374994
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
