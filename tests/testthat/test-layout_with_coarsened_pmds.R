se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE)

cg <- CellGraphs(se)[[1]]
g <- cg@cellgraph

test_that("layout_with_coarsened_pmds works as expected", {
  expect_no_error(xyz <- layout_with_coarsened_pmds(g, resolution = 0.5, n_iter = 3))

  expected_result <- structure(
    c(
      0.0301740141572202,
      -0.0207088590816486,
      1.03754672607807,
      1.24723542527862,
      -0.450040043961243,
      -0.496178932006319,
      0.571212032623079,
      0.510558661645138,
      0.297938718576735,
      0.27891625957249,
      0.329018037294113,
      0.343419920046724,
      -0.349764819249835,
      -0.445393984283404,
      -0.0818727010144745,
      -0.127602673448677,
      -0.378029597271313,
      -0.332940913185244
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
