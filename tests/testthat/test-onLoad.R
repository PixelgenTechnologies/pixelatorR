test_that("verbosity is set on package load", {
  expect_true(getOption("pixelatorR.verbose"))
})

test_that("arrow output directory is set on package load", {
  expect_true(dir.exists(getOption("pixelatorR.arrow_outdir")))
})

test_that(".onLoad works as expected", {
  expect_invisible(pixelatorR:::.onLoad())

  # Get current values
  arrow_outdir <- getOption("pixelatorR.arrow_outdir")
  pixelatorR_verbose <- getOption("pixelatorR.verbose")

  # Set verbosity to NULL
  options(pixelatorR.verbose = NULL)
  expect_equal(getOption("pixelatorR.verbose"), NULL)
  expect_invisible(pixelatorR:::.onLoad())
  expect_true(getOption("pixelatorR.verbose"))

  # Set arrow dir to NULL
  options(pixelatorR.arrow_outdir = NULL)
  expect_equal(getOption("pixelatorR.arrow_outdir"), NULL)
  expect_invisible(pixelatorR:::.onLoad())
  expect_equal(getOption("pixelatorR.arrow_outdir"), file.path(getwd(), "edgelists"))

  # Restore options
  options(pixelatorR.verbose = pixelatorR_verbose, pixelatorR.arrow_outdir = arrow_outdir)
})
