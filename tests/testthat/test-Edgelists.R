library(dplyr)

pxl_file <- minimal_pna_pxl_file()

tmp_dir <- file.path(fs::path_temp(), "test")
fs::dir_create(tmp_dir)
for (i in 1:2) {
  fs::file_copy(
    minimal_pna_pxl_file(),
    file.path(tmp_dir, glue::glue("S{i}.pxl")),
    overwrite = TRUE
  )
}
pxl_files <- list.files(tmp_dir, full.names = TRUE)

# Create a merged Seurat object
se_list <- lapply(seq_along(pxl_files), function(i) {
  se <- ReadPNA_Seurat(pxl_files[i])
  se$sample_id <- glue::glue("Sample{i}")
  return(se)
})
merged_seurat_obj <- merge(se_list[[1]], se_list[-1], add.cell.ids = LETTERS[1:2])

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  test_that("Edgelists method works as expected", {

    # TODO: handle the case where the test is run on Windows
    skip_on_os("windows")
    expect_no_error(el <- Edgelists(merged_seurat_obj, lazy = FALSE))
    expect_true(inherits(el, "tbl_df"))
    expect_equal(dim(el), c(1057188, 7))

    # NOTE: lazy loading can be problematic in tests because we're opening a connection
    # to the same PXL file(s) multiple times. This should not be a huge issue in practice
    # because one would normally not open multiple connections to the same file. This
    # doesn't seem to be an issue in Mac OS and Linux, but it is on Windows which will
    # complain if there are multiple connections to the same file. Since the connections
    # are closed when the garbage collector runs, it is quite difficult to control when
    # the connections are closed. Here we explicitly close the connection to avoid the
    # error which would typically not be done in practice.
    expect_no_error(el <- Edgelists(merged_seurat_obj, meta_data_columns = "sample_id"))
    expect_true(inherits(el, "tbl_lazy"))
    expect_true("sample_id" %in% colnames(el))
    DBI::dbDisconnect(el$src$con)
  })

  test_that("Edgelists method fails with invalid input", {
    expect_error(el <- Edgelists("Invalid"))
    expect_error(el <- Edgelists(merged_seurat_obj, lazy = "Invalid"))
    expect_error(el <- Edgelists(merged_seurat_obj, meta_data_columns = "Invalid"))
    expect_error(el <- Edgelists(merged_seurat_obj, assay = "Invalid"))
  })
}
