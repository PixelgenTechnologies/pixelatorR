se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE)

cg <- CellGraphs(se)[[1]]
g <- cg@cellgraph

test_that("layout_with_coarsened_pmds works as expected", {
  expect_no_error(xyz <- layout_with_coarsened_pmds(g, resolution = 0.5, n_iter = 3))

  expected_result <- structure(
    c(
      -0.15860076495787,
      0.830664133601459,
      -0.749051461819469,
      0.39850907378186,
      0.186144764689264,
      1.06582135287169,
      -0.591810526257705,
      -0.388041415546794,
      -0.3001100102262,
      -0.647841130988018,
      0.216512616002797,
      -0.151772710941226,
      -0.301869770556259,
      -0.0403459530893663,
      -0.35209224761166,
      0.208119803571115,
      0.445261628564128,
      -0.0110564590011096
    ),
    dim = c(6L, 3L),
    dimnames = list(NULL, c("x", "y", "z"))
  )

  expect_equal(xyz %>% head(), expected_result, tolerance = 1e-6)
})

test_that("layout_with_coarsened_pmds fails with invalid input", {
  expect_error(layout_with_coarsened_pmds("Invalid"))
  expect_error(layout_with_coarsened_pmds(g, dim = 4))
  expect_error(layout_with_coarsened_pmds(g, dim = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, resolution = 0))
  expect_error(layout_with_coarsened_pmds(g, resolution = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, pivots = 5))
  expect_error(layout_with_coarsened_pmds(g, pivots = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, n_iter = -1))
  expect_error(layout_with_coarsened_pmds(g, n_iter = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, jitter_sd = 1))
  expect_error(layout_with_coarsened_pmds(g, jitter_sd = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, weight_edges_by = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, leiden_iterations = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, leiden_weighted = "Invalid"))
})
