pxl_file <- minimal_pna_pxl_file()

test_that("ReadPNA_proximity works as expected", {
  # tbl_df
  expect_no_error(proximity <- ReadPNA_proximity(pxl_file, verbose = FALSE))
  expect_identical(dim(proximity), c(58696L, 9L))
  expected_names <- c(
    "marker_1", "marker_2", "join_count", "join_count_expected_mean",
    "join_count_expected_sd", "join_count_z", "join_count_p", "component",
    "log2_ratio"
  )
  expect_equal(expected_names, names(proximity))
  expect_s3_class(proximity, "tbl_df")

  # tbl_lazy
  expect_no_error(proximity <- suppressWarnings(ReadPNA_proximity(pxl_file, lazy = TRUE, verbose = FALSE)))
  expect_s3_class(proximity, "tbl_lazy")
  expect_equal(expected_names, colnames(proximity))
})

test_that("ReadPNA_proximity fails with invalid input", {
  expect_error(proximity <- ReadPNA_proximity("Invalid"))
  expect_error(proximity <- ReadPNA_proximity(pxl_file, verbose = "Invalid"))
  expect_error(proximity <- ReadPNA_proximity(pxl_file, calc_log2ratio = "Invalid"))
  expect_error(proximity <- ReadPNA_proximity(pxl_file, return_tibble = "Invalid"))
  expect_error(proximity <- ReadPNA_proximity(pxl_file, lazy = "Invalid"))
})
