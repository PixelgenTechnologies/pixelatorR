se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], verbose = FALSE)

cg <- CellGraphs(se)[[1]]
g <- cg@cellgraph

test_that("layout_with_spectral works as expected", {
  skip_if_not_installed("irlba")
  skip_if_not_installed("RSpectra")

  # Default SVD path
  expect_no_error(xyz <- layout_with_spectral(g, seed = 123))

  expected_result <- structure(
    c(
      0.166130349672414,
      -0.666892358787575,
      0.565997417190045,
      -0.39462441303961,
      -0.349597378534834,
      -0.881565222521497,
      -0.875178048777611,
      -0.421029602284673,
      -0.40260258513163,
      -0.424615779864241,
      0.433291489168997,
      -0.132562070557414,
      -0.251970582532076,
      -0.117536409783481,
      0.146760246388789,
      -0.44202702205321,
      -0.235340040191232,
      0.0844719663422277
    ),
    dim = c(6L, 3L),
    dimnames = list(
      c(
        "61208583141770358",
        "69733109123764664",
        "16002757515879905",
        "59822389138925142",
        "64270251753030037",
        "11111585952318758"
      ),
      c("x", "y", "z")
    )
  )

  # irlba's partial SVD moves by about 1e-4 between releases (2.3.7 vs 2.4.1)
  expect_equal(abs(xyz %>% head()), abs(expected_result), tolerance = 1e-4)
  expect_equal(nrow(xyz), length(g))
  expect_equal(colnames(xyz), c("x", "y", "z"))

  expect_no_error(xyz2 <- layout_with_spectral(g, dim = 2, seed = 123))
  expect_equal(ncol(xyz2), 2L)
  expect_equal(colnames(xyz2), c("x", "y"))

  # Eigen unnormalized path
  expect_no_error(xyz <- layout_with_spectral(g, normalize_laplacian = FALSE, solver = "eigen", seed = 123))

  expected_result <- structure(
    c(
      -0.171358949791906,
      0.6720640910113,
      -0.593827152715825,
      0.407084025576413,
      0.371802425000281,
      0.880701472922028,
      -0.760246551414652,
      -0.398790651683933,
      -0.404326498642815,
      -0.310689885796998,
      0.456837359713467,
      -0.211405088929061,
      -0.434784286797812,
      -0.162770480447131,
      0.0285781686542876,
      -0.478108651566909,
      -0.116705503235333,
      0.0820329916544583
    ),
    dim = c(6L, 3L),
    dimnames = list(NULL, c("x", "y", "z"))
  )

  expect_equal(abs(xyz %>% head()), abs(expected_result), tolerance = 1e-6)

  # Eigen normalized path
  expect_no_error(xyz <- layout_with_spectral(g, solver = "eigen", seed = 123))

  expected_result <- structure(
    c(
      0.166136074566168,
      -0.666915340083127,
      0.566016921611108,
      -0.394638011905492,
      -0.349609425757878,
      -0.881595601495921,
      0.875208162880329,
      0.421043552759889,
      0.402615693630875,
      0.424630340282078,
      -0.433307543481484,
      0.132568092094278,
      0.251971568783767,
      0.117482410048834,
      -0.146840643092578,
      0.442029873680554,
      0.235230758731725,
      -0.0843251762114535
    ),
    dim = c(6L, 3L),
    dimnames = list(NULL, c("x", "y", "z"))
  )

  expect_equal(abs(xyz %>% head()), abs(expected_result), tolerance = 1e-6)
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
