options(pixelatorR.arrow_outdir = tempdir())

test_that("Data loading with ReadMPX_Seurat works", {
  pg_data <-
    ReadMPX_Seurat(
      system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"), overwrite = TRUE
    )
  expect_s4_class(pg_data, "Seurat")
  expect_equal(names(pg_data@assays), c("mpxCells"))
  pg_data <-
    ReadMPX_Seurat(
      system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
      overwrite = TRUE, load_cell_graphs = TRUE
    )
  expect_s4_class(pg_data, "Seurat")
})

test_that("Data loading fails when an invalid file format is provided", {
  expect_error(
    ReadMPX_Seurat(
      "Sample01_test.pixl"
    ),
    regexp = "doesn't exist"
  )
})
