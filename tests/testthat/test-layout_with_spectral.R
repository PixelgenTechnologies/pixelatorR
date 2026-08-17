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
      0.185402329261988,
      0.666892358787576,
      0.684281597132036,
      0.565997417190045,
      0.575268497598157,
      0.875178048777611,
      0.900113153701813,
      0.421029602284673,
      0.431991710578858,
      0.40260258513163,
      0.405222438538549,
      0.251970582532076,
      0.167599110925018,
      0.117536409783481,
      0.109001478266206,
      0.146760246388789,
      0.0460186862773287
    ),
    dim = c(6L, 3L),
    dimnames = list(
      c(
        "61208583141770358-umi1",
        "50526950249468550-umi2",
        "69733109123764664-umi1",
        "43235234960499656-umi2",
        "16002757515879905-umi1",
        "4606209975865882-umi2"
      ),
      c("x", "y", "z")
    )
  )

  expect_equal(abs(xyz %>% head()), abs(expected_result), tolerance = 1e-6)
  expect_equal(nrow(xyz), length(g))
  expect_equal(colnames(xyz), c("x", "y", "z"))

  expect_no_error(xyz2 <- layout_with_spectral(g, dim = 2, seed = 123))
  expect_equal(ncol(xyz2), 2L)
  expect_equal(colnames(xyz2), c("x", "y"))

  # Eigen unnormalized path
  expect_no_error(xyz <- layout_with_spectral(g, normalize_laplacian = FALSE, solver = "eigen", seed = 123))

  expected_result <- structure(
    c(
      0.171358949792271,
      0.193857132967981,
      0.672064091010292,
      0.688354279187052,
      0.593827152715712,
      0.593039596069279,
      0.760246551157201,
      0.807963805056796,
      0.398790651582335,
      0.412327966225161,
      0.404326498663115,
      0.375929801226946,
      0.434784287255888,
      0.36406845310238,
      0.162770480671218,
      0.156233300732135,
      0.0285781684027663,
      0.0705112977261724
    ),
    dim = c(6L, 3L),
    dimnames = list(NULL, c("x", "y", "z"))
  )

  expect_equal(abs(xyz %>% head()), abs(expected_result), tolerance = 1e-6)

  # Eigen normalized path
  expect_no_error(xyz <- layout_with_spectral(g, solver = "eigen", seed = 123))

  expected_result <- structure(
    c(
      0.166136074566108,
      0.185408718272943,
      0.666915340083099,
      0.68430517766539,
      0.566016921610942,
      0.575288321504512,
      0.87520816287973,
      0.900144030436831,
      0.421043552759598,
      0.432006073469487,
      0.402615693630528,
      0.405236104403502,
      0.251971568792407,
      0.167588195789632,
      0.117482410054247,
      0.108950842216591,
      0.146840643092429,
      0.0460492780105566
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
