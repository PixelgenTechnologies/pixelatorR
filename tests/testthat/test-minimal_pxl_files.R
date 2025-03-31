test_that("minimal_mpx_pxl_file works as expected", {
  expect_no_error(pxl_file <- minimal_mpx_pxl_file())
  expect_true(fs::file_exists(pxl_file))
  expect_true(fs::path_ext(pxl_file) == "pxl")
})

test_that("minimal_pna_pxl_file works as expected", {
  expect_no_error(pxl_file <- minimal_pna_pxl_file())
  expect_true(fs::file_exists(pxl_file))
  expect_true(fs::path_ext(pxl_file) == "pxl")
})
