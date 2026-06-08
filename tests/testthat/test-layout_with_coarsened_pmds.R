se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE)

cg <- CellGraphs(se)[[1]]
g <- cg@cellgraph

test_that("layout_with_coarsened_pmds works as expected", {
  expect_no_error(xyz <- layout_with_coarsened_pmds(g, resolution = 0.1))

  expected_result <- structure(
    c(
      0.185070167476554,
      0.215214478324424,
      -0.977301209597777,
      -1.10128628386407,
      0.894203207119398,
      0.894454605209395,
      0.445430754370736,
      0.398746812263141,
      0.108245062416732,
      0.132919733209312,
      0.196100576528908,
      0.208913709744924,
      0.0520992694411212,
      -0.0372731279687765,
      0.0716099953647922,
      -0.0391260478142113,
      -0.338452488736382,
      -0.196268460201595
    ),
    dim = c(6L, 3L),
    dimnames = list(NULL, c("x", "y", "z"))
  )

  expect_equal(xyz %>% head(), expected_result, tolerance = 1e-6)

  expect_no_error(xyz <- layout_with_coarsened_pmds(g, resolution = 0.1, weight_edges_by = "crossing_edges"))

  expected_result <- structure(
    c(
      -0.220274499289863,
      -0.164322788816015,
      -1.06255147661951,
      -1.10104770581197,
      0.473308708048199,
      0.464284026496892,
      -0.0133807725077418,
      -0.0869502068496076,
      -0.0676139496428066,
      -0.073140660718092,
      -0.00579703502445285,
      0.0349425190088179,
      -0.507355237915104,
      -0.442826710702652,
      0.0987981981416092,
      0.190674453275647,
      -0.662302994281157,
      -0.548421628378492
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
