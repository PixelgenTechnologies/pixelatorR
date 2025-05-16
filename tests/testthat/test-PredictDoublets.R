
pxl_file <- minimal_pna_pxl_file()

seur <- ReadPNA_Seurat(pxl_file, verbose = FALSE)

test_that("PredictDoublets works as expected", {
  expect_no_error(PredictDoublets(seur))
  expect_no_error(suppressWarnings(PredictDoublets(seur, simulation_rate = 100)))

  # Simulate data
  set.seed(37)
  sim_data <-
    rbind(
      c(rnorm(n = 1000, mean = 0, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 0, sd = 1), rnorm(n = 1000, mean = 3, sd = 1)),
      c(rnorm(n = 1000, mean = 2, sd = 1), rnorm(n = 1000, mean = 2, sd = 1)),
      c(rnorm(n = 1000, mean = 4, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 3, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 2, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 0, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 0, sd = 1), rnorm(n = 1000, mean = 3, sd = 1)),
      c(rnorm(n = 1000, mean = 2, sd = 1), rnorm(n = 1000, mean = 2, sd = 1)),
      c(rnorm(n = 1000, mean = 4, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 3, sd = 1), rnorm(n = 1000, mean = 4, sd = 1)),
      c(rnorm(n = 1000, mean = 2, sd = 1), rnorm(n = 1000, mean = 4, sd = 1))
    ) %>%
    exp()
  sim_data <- as(sim_data, "CsparseMatrix")

  rownames(sim_data) <- paste0("CD", seq_len(nrow(sim_data)))
  colnames(sim_data) <- paste0("cell", seq_len(ncol(sim_data)))

  expect_no_error(pred_seur <- PredictDoublets(sim_data))

  expect_named(
    pred_seur,
    c(
      "doublet_nns", "doublet_nn_rate",
      "doublet_p", "doublet_p_adj",
      "doublet_prediction"
    )
  )
  expect_true(all(paste0("cell", seq_len(ncol(sim_data))) %in% rownames(pred_seur)))
  expect_equal(
    head(pred_seur),
    structure(
      list(
        doublet_nns = c(
          76L,
          67L,
          56L,
          60L,
          84L,
          92L
        ),
        doublet_nn_rate = c(
          0.76,
          0.67,
          0.56,
          0.6,
          0.84,
          0.92
        ),
        doublet_p = c(
          0.461671132081412,
          0.972405435869905,
          0.999989078555761,
          0.999676034583683,
          0.0211106216250894,
          1.20854259813212e-05
        ),
        doublet_p_adj = c(
          0.842465569491628,
          0.999999999441245,
          0.999999999441245,
          0.999999999441245,
          0.0657651764021476,
          0.000199759107129277
        ),
        doublet_prediction = c(
          "singlet",
          "singlet",
          "singlet",
          "singlet",
          "singlet",
          "doublet"
        )
      ),
      row.names = c(
        "cell1",
        "cell10",
        "cell100",
        "cell1000",
        "cell1001",
        "cell1002"
      ),
      class = "data.frame"
    )
  )


  # Errors
  expect_error(PredictDoublets(sim_data, simulation_rate = -1))
  expect_error(PredictDoublets(sim_data, n_neighbor = -1))
  expect_error(PredictDoublets(sim_data, npcs = -1))
  expect_error(PredictDoublets(sim_data, p_adjust_method = "notamethod"))
  expect_error(PredictDoublets(sim_data, p_threshold = "1"))
  expect_error(PredictDoublets(sim_data, seed = "1"))
})
