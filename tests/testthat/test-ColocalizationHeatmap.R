prox <- ReadPNA_proximity(minimal_pna_pxl_file())
prox_summarized <- prox %>%
  slice_sample(n = 1e4) %>%
  group_by(marker_1, marker_2) %>%
  summarize(mean_log2_ratio = mean(log2_ratio), .groups = "drop") %>%
  mutate(test1 = marker_1, test2 = marker_2, estimate = mean_log2_ratio) %>%
  ungroup()

test_that("ColocalizationHeatmap works as expected", {
  # Default method
  expect_no_error(ColocalizationHeatmap(prox_summarized))

  # Dots method
  expect_no_error(p <- ColocalizationHeatmap(prox_summarized, type = "dots", size_col = "mean_log2_ratio"))
  expect_true(is.factor(p$data$marker_1))
  expect_true(is.factor(p$data$marker_2))
  expect_equal(p$data$marker_1 %>% levels(), p$data$marker_2 %>% levels())

  # Custom marker_1 and marker_2 columns
  expect_no_error(ColocalizationHeatmap(prox_summarized, marker1_col = "test1", marker2_col = "test2"))

  # Custom value column
  expect_no_error(ColocalizationHeatmap(prox_summarized, value_col = "mean_log2_ratio"))

  # Symmetrise FALSE
  expect_no_error(p <- ColocalizationHeatmap(prox_summarized, symmetrise = FALSE, type = "dots", size_col = "mean_log2_ratio"))
  expect_true(is.factor(p$data$marker_1))
  expect_true(is.factor(p$data$marker_2))
  expect_true(!identical(p$data$marker_1 %>% levels(), p$data$marker_2 %>% levels()))
})

test_that("ColocalizationHeatmap fails with invalid input", {
  expect_error(ColocalizationHeatmap("Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, marker1_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, marker2_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, value_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, type = "dots", size_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, size_col_transform = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, size_range = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, colors = FALSE))
  expect_error(ColocalizationHeatmap(prox_summarized, cluster_rows = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, cluster_cols = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, type = "Invalid"))
})
