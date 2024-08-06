pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file) %>%
  LoadCellGraphs(cells = colnames(.)[2])
g <- CellGraphs(seur_obj)[[2]] %>%
  CellGraphData("cellgraph")


test_that("layout_with_weighted_pmds works as expected", {
  set.seed(123)

  # Compute weighted layout (cos_dist)
  expect_no_error({
    layout <- g %>%
      layout_with_weighted_pmds(pivots = 200, dim = 3, method = "cos_dist") %>%
      as_tibble(.name_repair = ~ c("x", "y", "z"))
  })
  layout_expected <- structure(
    list(
      x = c(
        518.739441082951,
        -304.40481679007,
        486.359827414132, -449.579623921567,
        953.163978220311,
        -788.780148291502
      ),
      y = c(
        -172.069727647667, -360.954231392934,
        -383.634163391642,
        386.703832810495,
        73.1979429546171,
        42.748662559004
      ),
      z = c(
        362.12952104792,
        -380.060102715779,
        348.334211929446,
        421.718168743799,
        -157.936334425518,
        -41.2389861032635
      )
    ),
    row.names = c(NA, -6L),
    class = c("tbl_df", "tbl", "data.frame")
  )
  expect_equal(head(layout) %>% mutate(across(everything(), ~ abs(.x))), layout_expected %>% mutate(across(everything(), ~ abs(.x))))

  # Compute weighted layout (prob_dist)
  expect_no_error({
    layout <- g %>%
      layout_with_weighted_pmds(pivots = 200, dim = 3, method = "prob_dist") %>%
      as_tibble(.name_repair = ~ c("x", "y", "z"))
  })
  layout_expected <- structure(
    list(
      x = c(
        -13826790.7416902,
        3518440.08631891,
        4454564.71551799,
        9384512.60045416,
        -17927509.7827884,
        27080144.2133983
      ),
      y = c(
        19958497.8955713,
        19055.4101331546,
        25416819.1467022,
        -4171319.32211823,
        -10872496.9390494, -14700709.8559556
      ),
      z = c(
        -7479934.55627813,
        16819433.9054322, -1105173.82157912,
        -21954674.4342814,
        4014742.5689878,
        -4560062.21766825
      )
    ),
    row.names = c(NA, -6L),
    class = c("tbl_df", "tbl", "data.frame")
  )
  expect_equal(head(layout) %>% mutate(across(everything(), ~ abs(.x))), layout_expected %>% mutate(across(everything(), ~ abs(.x))))
})

test_that("layout_with_weighted_pmds fails with invalid input", {
  # 0 pivots
  expect_error(layout <- layout_with_weighted_pmds(g, pivots = 0))

  # Invalid dim
  expect_error(layout <- layout_with_weighted_pmds(g, pivots = 200, dim = 1))

  # Invalid method
  expect_error(layout <- layout_with_weighted_pmds(g, pivots = 200, dim = 3, method = "Invalid"))
})
