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

  TEMP <- Sys.getenv("TMPDIR")
  session_tmp <- tempdir()
  # Temporarily change TMPDIR
  temp_dir_use <- file.path(TEMP, "pixelatorR_test_tmp")
  dir.create(temp_dir_use, showWarnings = FALSE)
  Sys.setenv(TMPDIR = temp_dir_use)
  unlink(tempdir(), recursive = TRUE)
  cur_temp_dir <- tempdir(check = TRUE)

  # Tempdir should be empty
  expect_no_error({
    pg_data <-
      ReadMPX_counts(
        system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")
      )
  })
  expect_true(length(list.files(cur_temp_dir)) == 0)
  Sys.setenv(TMPDIR = TEMP)
  unlink(temp_dir_use, recursive = TRUE)
  cur_temp_dir <- tempdir(check = TRUE)
})
