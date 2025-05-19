pxl_file <- minimal_pna_pxl_file()

seur <- ReadPNA_Seurat(pxl_file, verbose = FALSE)

test_that("PredictDoublets works as expected", {
  expect_no_error(PredictDoublets(seur))
  expect_no_error(suppressWarnings(PredictDoublets(seur, simulation_rate = 100)))

  # Simulate data
  set.seed(37)
  sim_data <-
    rbind(
      c(rnorm(n = 200, mean = 0, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 0, sd = 1), rnorm(n = 200, mean = 3, sd = 1)),
      c(rnorm(n = 200, mean = 2, sd = 1), rnorm(n = 200, mean = 2, sd = 1)),
      c(rnorm(n = 200, mean = 4, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 3, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 2, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 0, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 0, sd = 1), rnorm(n = 200, mean = 3, sd = 1)),
      c(rnorm(n = 200, mean = 2, sd = 1), rnorm(n = 200, mean = 2, sd = 1)),
      c(rnorm(n = 200, mean = 4, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 3, sd = 1), rnorm(n = 200, mean = 4, sd = 1)),
      c(rnorm(n = 200, mean = 2, sd = 1), rnorm(n = 200, mean = 4, sd = 1))
    ) %>%
    exp()
  sim_data <- as(sim_data, "CsparseMatrix")

  rownames(sim_data) <- paste0("CD", seq_len(nrow(sim_data)))
  colnames(sim_data) <- paste0("cell", seq_len(ncol(sim_data)))

  expect_no_error(pred_seur <- PredictDoublets(sim_data, n_neighbor = 10))

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
        doublet_nns = c(7L, 5L, 9L, 7L, 8L, 8L),
        doublet_nn_rate = c(0.7, 0.5, 0.9, 0.7, 0.8, 0.8),
        doublet_p = c(
          0.775875091552734, 0.98027229309082,
          0.244025230407715, 0.775875091552734,
          0.525592803955078, 0.525592803955078
        ),
        doublet_p_adj = c(
          0.910117409446023, 0.995200297554132,
          0.53048963132112, 0.910117409446023,
          0.748174809900467, 0.748174809900467
        ),
        doublet_prediction = c(
          "singlet", "singlet", "singlet",
          "singlet", "singlet", "singlet"
        )
      ),
      row.names = c(
        "cell1", "cell2", "cell3",
        "cell4", "cell5", "cell6"
      ),
      class = "data.frame"
    )
  )


  expect_no_error(pred_seur <- PredictDoublets(sim_data,
    simulation_rate = 0.1,
    n_neighbor = 10,
    ref_cells1 = 1:10,
    ref_cells2 = 201:220
  ))
  expect_no_error(pred_seur <- PredictDoublets(sim_data,
    simulation_rate = 0.1,
    n_neighbor = 10,
    ref_cells1 = colnames(sim_data)[1:10],
    ref_cells2 = colnames(sim_data)[201:220]
  ))


  # Errors
  expect_error(PredictDoublets(sim_data, simulation_rate = -1))
  expect_error(PredictDoublets(sim_data, n_neighbor = -1))
  expect_error(PredictDoublets(sim_data, npcs = -1))
  expect_error(PredictDoublets(sim_data, p_adjust_method = "notamethod"))
  expect_error(PredictDoublets(sim_data, p_threshold = "1"))
  expect_error(PredictDoublets(sim_data, seed = "1"))
})
