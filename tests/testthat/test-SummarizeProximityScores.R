pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file, verbose = FALSE)
seur_obj$sample_id <- c("A", "A", "B", "B", "B")
proximity <- ProximityScores(
  seur_obj,
  add_marker_counts = TRUE,
  add_marker_proportions = TRUE,
  meta_data_columns = "sample_id"
)
proximity_lazy <- ProximityScores(
  seur_obj,
  add_marker_counts = TRUE,
  add_marker_proportions = TRUE,
  lazy = TRUE,
  meta_data_columns = "sample_id"
)

test_that("SummarizeProximityScores works as expected", {
  expect_no_error(proximity_summarized <- SummarizeProximityScores(proximity))
  expect_no_error(
    proximity_summarized_lazy <- SummarizeProximityScores(proximity_lazy)
  )
  expect_equal(dim(proximity_summarized), dim(proximity_summarized_lazy))

  # Using log2_ratio
  expect_no_error(
    proximity_summarized <- SummarizeProximityScores(
      proximity,
      group_vars = "sample_id"
    )
  )
  expect_no_error(
    proximity_summarized_lazy <- SummarizeProximityScores(
      proximity_lazy,
      group_vars = "sample_id"
    )
  )
  expect_equal(dim(proximity_summarized), dim(proximity_summarized_lazy))
  expect_equal(dim(proximity_summarized), c(25122, 8))

  # Using join_count_z
  expect_no_error(
    proximity_summarized <- SummarizeProximityScores(
      proximity,
      proximity_metric = "join_count_z"
    )
  )
  expect_no_error(
    proximity_summarized_lazy <- SummarizeProximityScores(
      proximity_lazy,
      proximity_metric = "join_count_z"
    )
  )
  expect_equal(dim(proximity_summarized), dim(proximity_summarized_lazy))
  expect_equal(dim(proximity_summarized), c(12561, 7))

  # Including or excluding missing obs
  expect_no_error(
    proximity_summarized_w0 <- SummarizeProximityScores(
      proximity,
      include_missing_obs = TRUE,
      detailed = TRUE
    )
  )
  expect_no_error(
    proximity_summarized_wo0 <- SummarizeProximityScores(
      proximity,
      include_missing_obs = FALSE,
      detailed = TRUE
    )
  )
  expect_true(all(sapply(proximity_summarized_w0$log2_ratio_list, length) == 5))
  expect_true(
    !all(sapply(proximity_summarized_wo0$log2_ratio_list, length) == 5)
  )
  expect_true(
    !all(
      proximity_summarized_w0$mean_log2_ratio ==
        proximity_summarized_wo0$mean_log2_ratio
    )
  )

  # Non-detailed mean and median match detailed list calculations
  for (inc_miss in c(TRUE, FALSE)) {
    for (stat in c("mean", "median")) {
      res_simple <- SummarizeProximityScores(
        proximity,
        group_vars = "sample_id",
        include_missing_obs = inc_miss,
        summary_stat = stat,
        detailed = FALSE
      ) |>
        arrange(sample_id, marker_1, marker_2)
      res_det <- SummarizeProximityScores(
        proximity,
        group_vars = "sample_id",
        include_missing_obs = inc_miss,
        summary_stat = stat,
        detailed = TRUE
      ) |>
        arrange(sample_id, marker_1, marker_2)
      col_name <- glue::glue("{stat}_log2_ratio")
      expect_equal(res_simple[[col_name]], res_det[[col_name]])
      expect_equal(
        colnames(res_simple),
        c(
          "sample_id",
          "marker_1",
          "marker_2",
          "n_cells_detected",
          "n_cells",
          "n_cells_missing",
          "pct_detected",
          col_name
        )
      )
    }
  }

  # Multiple proximity metrics at once
  res_both <- SummarizeProximityScores(
    proximity,
    proximity_metric = c("log2_ratio", "join_count_z"),
    group_vars = "sample_id"
  ) |>
    arrange(sample_id, marker_1, marker_2)
  res_log2 <- SummarizeProximityScores(
    proximity,
    proximity_metric = "log2_ratio",
    group_vars = "sample_id"
  ) |>
    arrange(sample_id, marker_1, marker_2)
  res_jcz <- SummarizeProximityScores(
    proximity,
    proximity_metric = "join_count_z",
    group_vars = "sample_id"
  ) |>
    arrange(sample_id, marker_1, marker_2)
  expect_equal(res_both$mean_log2_ratio, res_log2$mean_log2_ratio)
  expect_equal(res_both$mean_join_count_z, res_jcz$mean_join_count_z)
})

test_that("SummarizeProximityScores propagates observed NAs and matches NA group_vars", {
  prox_edge <- tibble::tibble(
    component = c("c1", "c2", "c3", "c1", "c2", "c4", "c5"),
    sample_id = c("A", "A", "A", "A", "A", NA_character_, NA_character_),
    marker_1 = c("M1", "M1", "M2", "M2", "M2", "M1", "M2"),
    marker_2 = c("M1", "M1", "M2", "M2", "M2", "M1", "M2"),
    join_count = c(1, 2, 1, 2, 3, 2, 1),
    join_count_expected_mean = c(1, 1, 1, 1, 1, 1, 1),
    log2_ratio = c(2.0, NA_real_, 1.0, 2.0, 4.0, 3.0, NA_real_),
    join_count_z = c(1.5, NA_real_, 0.5, 1.0, 2.0, 4.0, 2.0)
  )

  for (inc_miss in c(TRUE, FALSE)) {
    for (stat in c("mean", "median")) {
      res_simple <- SummarizeProximityScores(
        prox_edge,
        proximity_metric = c("log2_ratio", "join_count_z"),
        group_vars = "sample_id",
        include_missing_obs = inc_miss,
        summary_stat = stat,
        detailed = FALSE
      ) |>
        arrange(sample_id, marker_1, marker_2)

      res_det <- SummarizeProximityScores(
        prox_edge,
        proximity_metric = c("log2_ratio", "join_count_z"),
        group_vars = "sample_id",
        include_missing_obs = inc_miss,
        summary_stat = stat,
        detailed = TRUE
      ) |>
        arrange(sample_id, marker_1, marker_2)

      col_log2 <- glue::glue("{stat}_log2_ratio")
      col_jcz <- glue::glue("{stat}_join_count_z")

      expect_equal(res_simple$n_cells, res_det$n_cells)
      expect_equal(res_simple$n_cells_missing, res_det$n_cells_missing)
      expect_equal(res_simple$pct_detected, res_det$pct_detected)
      expect_equal(res_simple[[col_log2]], res_det[[col_log2]])
      expect_equal(res_simple[[col_jcz]], res_det[[col_jcz]])

      # Group ("A", "M1", "M1") has c(2.0, NA_real_) and 1 missing cell -> must be NA
      row_a_m1 <- filter(res_simple, sample_id == "A", marker_1 == "M1")
      expect_equal(row_a_m1[[col_log2]], NA_real_)
      expect_equal(row_a_m1[[col_jcz]], NA_real_)

      # Group (NA, "M1", "M1") has 1 detected cell (3.0) out of 2 total NA-sample cells
      row_na_m1 <- filter(res_simple, is.na(sample_id), marker_1 == "M1")
      expect_equal(row_na_m1$n_cells, 2L)
      expect_equal(row_na_m1$n_cells_missing, 1L)
      expect_equal(row_na_m1$pct_detected, 0.5)
      expect_equal(row_na_m1[[col_log2]], if (inc_miss) 1.5 else 3.0)
    }
  }
})

test_that("SummarizeProximityScores fails with invalid input", {
  expect_error(proximity_filtered <- SummarizeProximityScores("Invalid"))
  expect_error(
    proximity_filtered <- SummarizeProximityScores(
      proximity,
      proximity_metric = "Invalid"
    )
  )
  expect_error(
    proximity_filtered <- SummarizeProximityScores(
      proximity,
      proximity_metric = "marker_1"
    )
  )
  expect_error(
    proximity_filtered <- SummarizeProximityScores(
      proximity,
      group_vars = "Invalid"
    )
  )
  expect_error(
    proximity_filtered <- SummarizeProximityScores(
      proximity,
      group_vars = "join_count"
    )
  )
  expect_error(
    proximity_filtered <- SummarizeProximityScores(bind_rows(
      proximity,
      proximity
    )),
    regexp = "(duplicate rows|multiple rows)"
  )
})
