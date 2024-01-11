el <- ReadMPX_arrow_edgelist(system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR"),
                       outdir = tempdir(),
                       overwrite = TRUE)

test_that("export_edgelist_to_parquet works as expected", {
  expect_no_error(export_edgelist_to_parquet(object = el, outdir = tempdir(), overwrite = TRUE))
})


test_that("export_edgelist_to_parquet fails when invalid input is provided", {
  expect_error(export_edgelist_to_parquet(object = "Invalid"))
  expect_error(export_edgelist_to_parquet(object = el, outdir = TRUE), "'outdir' must be a character vector of length 1")
  expect_error(suppressWarnings({export_edgelist_to_parquet(object = el, outdir = "__Invalid directory__")}), "outdir .* doesn't exist")
})
