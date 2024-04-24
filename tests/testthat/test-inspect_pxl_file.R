pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)

test_that("inspect_pxl_file works as expected", {

  expect_no_error({pxl_file_info <- inspect_pxl_file(pxl_file)})

  expect_identical(
    pxl_file_info,
    structure(
      list(
        file_type = c(
          "adata.h5ad",
          "edgelist.parquet",
          "metadata.json",
          "polarization.parquet",
          "colocalization.parquet"
        ),
        n = c(1L, 1L, 1L, 1L, 1L),
        file = list(
          "adata.h5ad",
          "edgelist.parquet",
          "metadata.json",
          "polarization.parquet",
          "colocalization.parquet"
        )
      ),
      class = c("tbl_df",
                "tbl", "data.frame"), row.names = c(NA, -5L)))
})

test_that("inspect_pxl_file fails with invalid input", {

  expect_error({pxl_file_info <- inspect_pxl_file("Invalid")}, "File 'Invalid' does not exist.")

})
