pxl_file <- minimal_pna_pxl_file()

test_that("ReadPNA_layouts works as expected", {
  # tbl_df
  expect_no_error(layouts <- ReadPNA_layouts(pxl_file, verbose = FALSE))
  expect_equal(length(layouts), 5)
  expected_names <- c("name", "x", "y", "z")
  expect_equal(expected_names, names(layouts[[1]]))
  expect_type(layouts, "list")
  expect_equal(dim(layouts[[1]]), c(43543, 4))
})

test_that("ReadPNA_proximity fails with invalid input", {
  expect_error(layouts <- ReadPNA_layouts("Invalid"))
  expect_error(layouts <- ReadPNA_layouts(pxl_file, verbose = "Invalid"))
  expect_error(layouts <- ReadPNA_layouts(pxl_file, cells = "Invalid"))
  expect_error(layouts <- ReadPNA_layouts(pxl_file, add_marker_counts = "Invalid"))
})
