se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE)

cg <- CellGraphs(se)[[1]]
g <- cg@cellgraph

test_that("layout_with_coarsened_pmds works as expected", {
  expect_no_error(xyz <- layout_with_coarsened_pmds(g, resolution = 0.5, n_iter = 3))

  expected_result <- structure(
    c(
      0.0175374496643101,
      -0.00995572874022552,
      1.02077048105682,
      1.07658885150354,
      -0.42351419724422,
      -0.470416679970342,
      0.54984973245313,
      0.517271714479005,
      0.312237935798172,
      0.310466441132375,
      0.343200818771736,
      0.328726989744251,
      -0.401195946972589,
      -0.455131491839664,
      -0.0855087613342495,
      -0.104659414138358,
      -0.359347957201346,
      -0.346398097162904
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
  expect_error(layout_with_coarsened_pmds(g, jitter_sd = 0))
  expect_error(layout_with_coarsened_pmds(g, jitter_sd = "Invalid"))
  expect_error(layout_with_coarsened_pmds(g, weight_edges_by = "Invalid"))
})
