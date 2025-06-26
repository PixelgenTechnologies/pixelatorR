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

  set.seed(1)
  expect_equal(SimulateDoublets(sim_data, n_sim = 2),
               new("dgCMatrix",
                   i = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
                         10L, 11L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L),
                   p = c(0L, 12L, 24L), Dim = c(12L, 2L),
                   Dimnames = list(c("CD1", "CD2", "CD3", "CD4", "CD5",
                                     "CD6", "CD7", "CD8", "CD9", "CD10",
                                     "CD11", "CD12"),
                                   c("simcell_1", "simcell_2")),
                   x = c(0.892377550966689,
                         1.9919700963907, 8.19334019906415, 33.9668761030975, 96.8857558686188,
                         6.73949417506845, 1.77609666447528, 5.72555905996148, 10.5020990941118,
                         196.929993085046, 48.4059579357175, 7.1179363685714, 13.2040746793764,
                         6.38029192043951, 10.9703459687392, 137.604049023378, 37.7697104966331,
                         194.054547743392, 34.7757627417018, 5.63197859017003, 4.18806413836322,
                         22.7274280516269, 97.7250280221664, 14.3866787100948), factors = list()))
  set.seed(1)
  expect_equal(SimulateDoublets(sim_data, n_sim = 2, return_id = TRUE),
               list(counts = new("dgCMatrix",
                                 i = c(0L, 1L, 2L, 3L, 4L, 5L,
                                       6L, 7L, 8L, 9L, 10L, 11L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                                       9L, 10L, 11L), p = c(0L, 12L, 24L), Dim = c(12L, 2L),
                                 Dimnames = list(
                                   c("CD1", "CD2", "CD3", "CD4", "CD5", "CD6", "CD7", "CD8",
                                     "CD9", "CD10", "CD11", "CD12"), c("simcell_1", "simcell_2"
                                     )),
                                 x = c(0.892377550966689, 1.9919700963907, 8.19334019906415,
                                       33.9668761030975, 96.8857558686188, 6.73949417506845, 1.77609666447528,
                                       5.72555905996148, 10.5020990941118, 196.929993085046, 48.4059579357175,
                                       7.1179363685714, 13.2040746793764, 6.38029192043951, 10.9703459687392,
                                       137.604049023378, 37.7697104966331, 194.054547743392, 34.7757627417018,
                                       5.63197859017003, 4.18806413836322, 22.7274280516269, 97.7250280221664,
                                       14.3866787100948),
                                 factors = list()),
                    ids = structure(list(i1 = c(182L,
                                                47L),
                                         i2 = c(24L, 306L),
                                         id1 = c("cell182",
                                                 "cell47"),
                                         id2 = c("cell24",
                                                 "cell306")),
                                    row.names = c(NA, -2L), class = c("tbl_df", "tbl",
                                                                      "data.frame"))))



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
