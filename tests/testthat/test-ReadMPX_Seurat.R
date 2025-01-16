options(Seurat.object.assay.version = "v3")

test_that("Data loading with ReadMPX_Seurat works", {
  pg_data <-
    ReadMPX_Seurat(
      system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
      overwrite = TRUE
    )
  expect_s4_class(pg_data, "Seurat")
  expect_equal(names(pg_data@assays), c("mpxCells"))
  pg_data <-
    ReadMPX_Seurat(
      system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
      overwrite = TRUE, load_cell_graphs = TRUE
    )
  expect_s4_class(pg_data, "Seurat")

  # Test that spatial scores are loaded correctly
  pg_data <-
    ReadMPX_Seurat(
      system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
      overwrite = TRUE,
      load_polarity_scores = FALSE,
      load_colocalization_scores = TRUE
    )

  expect_equal(nrow(pg_data@assays$mpxCells@polarization), 0)
  expect_equal(nrow(pg_data@assays$mpxCells@colocalization), 14649)

  pg_data <-
    ReadMPX_Seurat(
      system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR"),
      overwrite = TRUE,
      load_polarity_scores = TRUE,
      load_colocalization_scores = FALSE
    )

  expect_equal(nrow(pg_data@assays$mpxCells@polarization), 380)
  expect_equal(nrow(pg_data@assays$mpxCells@colocalization), 0)
})

test_that("Data loading fails when an invalid file format is provided", {
  expect_error(
    ReadMPX_Seurat(
      "Sample01_test.pixl"
    )
  )
})
