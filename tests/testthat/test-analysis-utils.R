plot <-
  ggplot()

test_that("export_plot works as expected", {
  temp_file <- tempfile()

  expect_no_error(
    export_plot(
      filename = temp_file,
      plot = plot
    )
  )

  expect_true(file.exists(paste0(temp_file, ".png")))
  expect_true(file.exists(paste0(temp_file, ".pdf")))

  temp_file <- tempfile()

  expect_no_error(
    export_plot(
      filename = paste0(temp_file, ".png"),
      plot = plot
    )
  )

  expect_true(file.exists(paste0(temp_file, ".png")))
  expect_false(file.exists(paste0(temp_file, ".pdf")))
})
