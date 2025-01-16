library(dplyr)

# Load example data as a Seurat object
pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)

# Copy PXL file to tempdir
tmp_pxl_file <- file.path(fs::path_temp(), "five_cells.pxl")
fs::file_copy(pxl_file, tmp_pxl_file)
seur_obj <- ReadMPX_Seurat(tmp_pxl_file)

# Now we can load graphs
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])

# Removing or moving PXL file will make graphs inaccessible
fs::file_delete(tmp_pxl_file)

test_that("RestorePaths works as expected", {
  pxl_files_dir <- system.file("extdata/five_cells",
    package = "pixelatorR"
  )
  expect_no_error(seur_obj <- RestorePaths(seur_obj, pxl_files_dir = pxl_files_dir))
  expect_no_error(seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1], force = TRUE))
})

test_that("RestorePaths fails with invalid input", {
  expect_error(seur_obj <- RestorePaths("Invalid"))
  expect_error(seur_obj <- RestorePaths(seur_obj, pxl_files_dir = "Invalid"))
})
