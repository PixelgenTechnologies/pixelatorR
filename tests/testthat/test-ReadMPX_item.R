testfile <- system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR")

test_that("ReadMPX_item works as expected", {

  # Read polarization scores
  polarization <- ReadMPX_item(filename = testfile, items = "polarization")
  expect_s3_class(polarization, "tbl_df")
  expect_equal(ncol(polarization), 6)
  expect_equal(nrow(polarization), 225)

  # Read colocalization scores
  colocalization <- ReadMPX_item(filename = testfile, items = "colocalization")
  expect_s3_class(colocalization, "tbl_df")
  expect_equal(ncol(colocalization), 15)
  expect_equal(nrow(colocalization), 2680)

  # Read edgelist
  edgelist <- ReadMPX_item(filename = testfile, items = "edgelist")
  expect_s3_class(edgelist, "tbl_df")
  expect_equal(ncol(edgelist), 9)
  expect_equal(nrow(edgelist), 64750)
})

test_that("ReadMPX_item fails when invalid input is provided", {
  expect_error(ReadMPX_item(filename = "invalid_file.pixl"), "Expected a .pxl file")
  expect_error(ReadMPX_item(filename = "invalid_file.pxl"), "doesn't exist")
  expect_error(ReadMPX_item(filename = testfile, items = "invalid_item"), "Invalid items")
})
