options(Seurat.object.assay.version = "v3")

test_that("Data loading with ReadMPX_counts works as expected", {
  expect_no_error({
    pg_data <-
      ReadMPX_counts(
        system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")
      )
  })
  expect_type(pg_data, "integer")
  expect_equal(dim(pg_data), c(80, 5))

  # Tempfile shouldn't be cleaned when return_list = TRUE
  expect_no_error({
    pg_data <-
      ReadMPX_counts(
        system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
        return_list = TRUE
      )
  })
  expect_type(pg_data$X, "integer")
  expect_equal(dim(pg_data$X), c(80, 5))
  expect_true(file.exists(pg_data$tmp_file))
  fs::file_delete(pg_data$tmp_file)
})
