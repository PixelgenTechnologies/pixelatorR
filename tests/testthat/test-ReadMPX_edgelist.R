test_that("ReadMPX_arrow_edgelist works as expected", {
  el <- ReadMPX_arrow_edgelist(system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR"),
                         outdir = tempdir(),
                         overwrite = TRUE)
  expect_type(el, type = "environment")
  expect_equal(names(el), c("upia","upib","marker","count",
                            "component","sample"))
  expect_no_error({msg <- capture_messages({el <- ReadMPX_arrow_edgelist(system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR"),
                         overwrite = TRUE)})})
})


test_that("ReadMPX_arrow_edgelist fails when invalid input is provided", {
  expect_error({el <- ReadMPX_arrow_edgelist("Invalid file",
                         outdir = tempdir())}, "file 'Invalid file' doesn't exist")
  expect_error({suppressWarnings({el <- ReadMPX_arrow_edgelist(system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR"),
                                       outdir = "__invalid directory__")})}, "outdir .* doesn't exist")
})
