pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file, verbose = FALSE)
seur_obj$sample_id <- c("A", "A", "B", "B", "B")
proximity <- ProximityScores(seur_obj, add_marker_counts = TRUE, add_marker_proportions = TRUE, meta_data_columns = "sample_id")
proximity_lazy <- ProximityScores(seur_obj, add_marker_counts = TRUE, add_marker_proportions = TRUE, lazy = TRUE, meta_data_columns = "sample_id")

test_that("SummarizeProximityScores works as expected", {
  expect_no_error(proximity_summarized <- SummarizeProximityScores(proximity))
  expect_no_error(proximity_summarized_lazy <- SummarizeProximityScores(proximity_lazy))
  expect_equal(dim(proximity_summarized), dim(proximity_summarized_lazy))

  # Using log2_ratio
  expect_no_error(proximity_summarized <- SummarizeProximityScores(proximity, group_vars = "sample_id"))
  expect_no_error(proximity_summarized_lazy <- SummarizeProximityScores(proximity_lazy, group_vars = "sample_id"))
  expect_equal(dim(proximity_summarized), dim(proximity_summarized_lazy))
  expect_equal(dim(proximity_summarized), c(25122, 8))

  # Using join_count_z
  expect_no_error(proximity_summarized <- SummarizeProximityScores(proximity, proximity_metric = "join_count_z"))
  expect_no_error(proximity_summarized_lazy <- SummarizeProximityScores(proximity_lazy, proximity_metric = "join_count_z"))
  expect_equal(dim(proximity_summarized), dim(proximity_summarized_lazy))
  expect_equal(dim(proximity_summarized), c(12561, 7))
})

test_that("SummarizeProximityScores fails with invalid input", {
  expect_error(proximity_filtered <- SummarizeProximityScores("Invalid"))
  expect_error(proximity_filtered <- SummarizeProximityScores(proximity, proximity_metric = "Invalid"))
  expect_error(proximity_filtered <- SummarizeProximityScores(proximity, proximity_metric = "marker_1"))
  expect_error(proximity_filtered <- SummarizeProximityScores(proximity, group_vars = "Invalid"))
  expect_error(proximity_filtered <- SummarizeProximityScores(proximity, group_vars = "join_count"))
})
