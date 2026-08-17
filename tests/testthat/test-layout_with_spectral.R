se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE)

cg <- CellGraphs(se)[[1]]
g <- cg@cellgraph

test_that("layout_with_spectral works with svd solver", {
  skip_if_not_installed("irlba")

  expect_no_error(xyz <- layout_with_spectral(g, solver = "svd", dim = 3, seed = 123))
  expect_equal(dim(xyz), c(igraph::vcount(g), 3))
  expect_equal(colnames(xyz), c("x", "y", "z"))
  expect_true(is.matrix(xyz))
  expect_true(all(is.finite(xyz)))

  expect_no_error(xyz2 <- layout_with_spectral(g, solver = "svd", dim = 2, seed = 123))
  expect_equal(dim(xyz2), c(igraph::vcount(g), 2))
  expect_equal(colnames(xyz2), c("x", "y"))
})

test_that("layout_with_spectral works with eigen solver", {
  skip_if_not_installed("RSpectra")

  expect_no_error(
    xyz <- layout_with_spectral(g, solver = "eigen", dim = 3, seed = 123)
  )
  expect_equal(dim(xyz), c(igraph::vcount(g), 3))
  expect_equal(colnames(xyz), c("x", "y", "z"))
  expect_true(all(is.finite(xyz)))

  expect_no_error(
    xyz_unnorm <- layout_with_spectral(
      g,
      solver = "eigen",
      normalize_laplacian = FALSE,
      dim = 3,
      seed = 123
    )
  )
  expect_equal(dim(xyz_unnorm), c(igraph::vcount(g), 3))
  expect_true(all(is.finite(xyz_unnorm)))
})

test_that("layout_with_spectral applies jitter when requested", {
  skip_if_not_installed("irlba")

  xyz <- layout_with_spectral(g, solver = "svd", dim = 3, seed = 123, jitter_sd = 0)
  xyz_jitter <- layout_with_spectral(g, solver = "svd", dim = 3, seed = 123, jitter_sd = 1e-2)
  expect_false(isTRUE(all.equal(xyz, xyz_jitter)))
})

test_that("layout_with_spectral fails with invalid input", {
  expect_error(layout_with_spectral("Invalid"))
  expect_error(layout_with_spectral(g, dim = 4))
  expect_error(layout_with_spectral(g, dim = "Invalid"))
  expect_error(layout_with_spectral(g, normalize_laplacian = "Invalid"))
  expect_error(layout_with_spectral(g, solver = "Invalid"))
  expect_error(layout_with_spectral(g, solver = "svd", normalize_laplacian = FALSE))
  expect_error(layout_with_spectral(g, jitter_sd = 1))
  expect_error(layout_with_spectral(g, jitter_sd = "Invalid"))
  expect_error(layout_with_spectral(g, seed = "Invalid"))
  expect_error(layout_with_spectral(g, verbose = "Invalid"))
})

test_that("ComputeLayout works with spectral layout method", {
  skip_if_not_installed("irlba")

  expect_no_error(layout <- ComputeLayout(g, layout_method = "spectral", dim = 3))
  expect_equal(nrow(layout), igraph::vcount(g))
  expect_equal(colnames(layout), c("x", "y", "z"))

  expect_no_error(cg_layout <- ComputeLayout(cg, layout_method = "spectral", dim = 3))
  expect_true("spectral_3d" %in% names(cg_layout@layout))
})
