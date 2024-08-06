seur <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
  overwrite = TRUE, return_cellgraphassay = TRUE
)

seur <- LoadCellGraphs(seur, cells = colnames(seur)[1])
seur <- seur %>% ComputeLayout(layout_method = "pmds", dim = 3, seed = 123)
cg <- CellGraphs(seur)[[1]]
g <- CellGraphData(cg, slot = "cellgraph")
counts <- CellGraphData(cg, slot = "counts")

test_that("local_G works as expected", {
  # Default settings
  expect_no_error(gi_mat <- local_G(g = g, counts = counts))
  expect_equal(dim(gi_mat), c(2470, 79))
  expect_type(gi_mat, "double")
  expected_values <- structure(
    c(
      -0.0918233864647522,
      5.60354339391612,
      4.97008915869393,
      1.60120673509987,
      0.12155828214979,
      -0.686429971341957,
      0.277856037783697,
      2.5315128050998,
      -0.118760619969203,
      -0.432928849830687,
      1.37658159346162, -0.164310530224204,
      -0.151069745965333,
      -0.2136795237295,
      -0.107213005054555,
      0.231100106817231,
      -0.651389682760557,
      1.45560324900578,
      0.0465557317169272, -0.335312458144784,
      -0.139016875779849,
      0.425590570987706,
      1.45472068367428, -0.325337007170695,
      -0.368374401123967
    ),
    dim = c(5L, 5L),
    dimnames = list(
      c(
        "GGTATATATTTTAAGTAGTTTAGTA-A",
        "ATTATGTTGTTGTATGTTTATTATT-A",
        "TCCGTGGTTTGATACTCGGGAATTT-A",
        "AGTGTAAGAGGTTGTTTCTTAGAAA-A",
        "GAGCAGACAATGGCGCTTAGCTAAA-A"
      ),
      c(
        "CD27", "HLA-ABC", "CD53",
        "B2M", "CD18"
      )
    )
  )
  expect_equal(gi_mat[1:5, 1:5], expected_values)

  # With p-values
  expect_no_error(gi_res <- local_G(g = g, counts = counts, return_p_vals = TRUE))
  expect_type(gi_res, "list")
  expect_equal(names(gi_res), c("gi_mat", "gi_p_mat", "gi_p_adj_mat"))
  expect_equal(gi_res$gi_mat[1:5, 1:5], expected_values)

  expected_p_values <-
    structure(
      c(
        0.926838362843666,
        2.10013507101039e-08,
        6.69221211840571e-07,
        0.109331137607092,
        0.903248854416789,
        0.49244201004635,
        0.78112287223606,
        0.0113571670485305,
        0.90546500852975,
        0.665066460234655,
        0.16864164491058,
        0.86948669243817,
        0.879920697785296,
        0.830796980251474,
        0.914619998357407,
        0.817237032136978,
        0.514794968962243,
        0.145502329370218,
        0.962867314717907,
        0.737389396939703,
        0.889436813591041,
        0.670406212142134,
        0.145746604049181,
        0.744926025272002,
        0.712594082456338
      ),
      dim = c(5L, 5L),
      dimnames = list(
        c(
          "GGTATATATTTTAAGTAGTTTAGTA-A",
          "ATTATGTTGTTGTATGTTTATTATT-A",
          "TCCGTGGTTTGATACTCGGGAATTT-A",
          "AGTGTAAGAGGTTGTTTCTTAGAAA-A",
          "GAGCAGACAATGGCGCTTAGCTAAA-A"
        ),
        c(
          "CD27", "HLA-ABC", "CD53",
          "B2M", "CD18"
        )
      )
    )

  expect_equal(gi_res$gi_p_mat[1:5, 1:5], expected_p_values)

  # Higher k
  expect_no_error(gi_mat <- local_G(g = g, counts = counts, k = 3))

  # k = 1 and use_weights=FALSE
  expect_no_error(gi_mat_use_weights_F <- local_G(g = g, counts = counts, use_weights = FALSE))

  # Higher k and use_weights=FALSE
  expect_no_error(gi_mat_use_k3_weights_F <- local_G(g = g, k = 3, counts = counts, use_weights = FALSE))

  # The two matrices should be different
  expect_true(!sum(gi_mat_use_weights_F == gi_mat_use_k3_weights_F))
})


test_that("local_G fails with invalid input expected", {
  # Invalid k
  expect_error(local_G(g = g, counts = counts, k = 0), "'k' must be and integer larger than 0")

  # Invalid counts
  expect_error(local_G(g = g, counts = "invalid"), "'counts' must be a sparse matrix of class 'dgCMatrix' or a 'matrix'")

  # Invalid g
  expect_error(local_G(g = "Invalid", counts = counts), "'g' must be an 'tbl_graph' or an 'igraph' object")

  expect_error(local_G(g = g, counts = counts, type = "Invalid"))
})
