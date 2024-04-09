test_that("verbosity is set on package load", {
  expect_true(getOption("pixelatorR.verbose"))
})

test_that(".onLoad works as expected", {
  expect_invisible(pixelatorR:::.onLoad())

  # Get current values
  pixelatorR_verbose <- getOption("pixelatorR.verbose")

  # Set verbosity to NULL
  options(pixelatorR.verbose = NULL)
  expect_equal(getOption("pixelatorR.verbose"), NULL)
  expect_invisible(pixelatorR:::.onLoad())
  expect_true(getOption("pixelatorR.verbose"))

  # Restore options
  options(pixelatorR.verbose = pixelatorR_verbose)
})
