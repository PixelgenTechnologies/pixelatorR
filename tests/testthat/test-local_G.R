seur <- ReadMPX_Seurat(system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
                       overwrite = TRUE, return_cellgraphassay = TRUE)

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
      -0.180818198677122,
      1.60447604867497,
      0.973660564365095,
      0.714869484955289,
      -0.0515512253415386,
      -1.15009352509223,
      -1.40111866040665,
      0.231963842286796,
      -0.89192459829829,
      1.96317708626894,
      2.17251720062701,-0.13428655212695,
      -0.123465217285206,
      -0.174634495200092,
      -0.0876222891636957,
      3.18250117577338,
      -1.01547397350661,
      -0.482733264613957,
      -0.218357325615071,-0.570503051696732,
      -0.380233943200656,
      0.205468985257245,
      -0.041929949219014,-0.449196887294857,
      -0.288389776392792
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
      c("CD27", "HLA-ABC", "CD53",
        "B2M", "CD18")
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
        0.856510278087473,
        0.108609159794192,
        0.330225104029602,
        0.474689694288205,
        0.958886284163075,
        0.250105352928739,
        0.161178592342013,
        0.816566097552408,
        0.37243334897577,
        0.0496256008740929,
        0.0298166736720564,
        0.893175986773189,
        0.901738716102041,
        0.861366841822041,
        0.930176886123393,
        0.00146008922961837,
        0.309879946800128,
        0.629285140857015,
        0.827150712482951,
        0.56833655289695,
        0.703771765019037,
        0.837205727318467,
        0.9665545413618,
        0.65328963273571,
        0.773048392717585
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
        c("CD27", "HLA-ABC", "CD53",
          "B2M", "CD18")
      )
    )

  expect_equal(gi_res$gi_p_mat[1:5, 1:5], expected_p_values)

  # Higher k
  expect_no_error(gi_mat <- local_G(g = g, counts = counts, k = 3))

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
