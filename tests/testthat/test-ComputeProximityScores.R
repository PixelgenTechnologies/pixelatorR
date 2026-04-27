library(dplyr)
se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
  LoadCellGraphs(cells = colnames(.)[1:2], verbose = FALSE)

cg_list <- CellGraphs(se)[1:2]
cg <- cg_list[[1]]

test_that("ComputeProximityScores works as expected", {
  # CellGraph
  expect_no_error(prox <- ComputeProximityScores(cg, mode = "permutation", iterations = 10))
  expect_no_error(ComputeProximityScores(cg, mode = "permutation", iterations = 10, calc_z_score = FALSE))
  expect_no_error(ComputeProximityScores(cg, mode = "permutation", iterations = 10, calc_log2_ratio = FALSE))

  expect_equal(
    structure(
      list(
        join_count = c(216, 8),
        join_count_expected_mean = c(68, 4.3),
        join_count_expected_sd = c(7.27247474309048, 0.948683298050514),
        marker_1 = c("B2M", "B2M"),
        marker_2 = c("B2M", "CD10"),
        join_count_z = c(20.3507066340263, 3.7),
        log2_ratio = c(1.66742466091313, 0.895663340185264)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    ),
    head(prox %>% arrange(marker_1, marker_2), 2)
  )

  old_locale <- Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  on.exit(Sys.setlocale("LC_COLLATE", old_locale), add = TRUE)
  expect_true(all(prox$marker_1 <= prox$marker_2))

  # Analytical
  expect_no_error(prox <- ComputeProximityScores(cg, mode = "analytical"))
  expect_equal(
    structure(
      list(
        join_count = c(216, 8),
        join_count_expected_mean = c(72.0819106992755, 4.10121216047602),
        join_count_expected_sd = c(6.83625825741441, 1.55421928931672),
        marker_1 = c("B2M", "B2M"),
        marker_2 = c("B2M", "CD10"),
        join_count_z = c(21.0521726771564, 2.50851849949565),
        log2_ratio = c(1.58332215361863, 0.963949622111638)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    ),
    head(prox %>% arrange(marker_1, marker_2), 2)
  )

  # list
  expect_no_error(ComputeProximityScores(cg_list, mode = "permutation", iterations = 10))
  # PNAAssay
  expect_no_error(ComputeProximityScores(se[["PNA"]], cells = colnames(se)[1:2], iterations = 10))

  # Seurat
  expect_no_error(ComputeProximityScores(se[["PNA"]], cells = colnames(se)[1:2], iterations = 10))
})

test_that("ComputeProximityScores fails with invalid input", {
  expect_error(ComputeProximityScores("Invalid"))
  expect_error(ComputeProximityScores(cg, iterations = "Invalid"))
  expect_error(ComputeProximityScores(cg, calc_z_score = "Invalid"))
  expect_error(ComputeProximityScores(cg, calc_log2_ratio = "Invalid"))
  expect_error(ComputeProximityScores(cg, seed = "Invalid"))
})
