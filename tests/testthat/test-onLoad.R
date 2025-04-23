test_that("verbosity is set on package load", {
  expect_true(getOption("pixelatorR.verbose"))
})

test_that(".onLoad works as expected", {
  # Get current values
  pixelatorR_verbose <- getOption("pixelatorR.verbose")

  # Set verbosity to NULL
  options(pixelatorR.verbose = NULL)
  expect_equal(getOption("pixelatorR.verbose"), NULL)

  # Restore options
  options(pixelatorR.verbose = pixelatorR_verbose)
})
