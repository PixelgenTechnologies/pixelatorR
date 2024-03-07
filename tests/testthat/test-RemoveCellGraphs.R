options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/mock_data",
                        "mock_mpx_data.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])

test_that("RemoveCellGraphs works for CellGraphAssay objects", {
  expect_no_error({cg_assay <- RemoveCellGraphs(seur_obj[["mpxCells"]])})
  classes <- sapply(cg_assay@cellgraphs, class)
  expect_true(all(classes %in% "NULL"))
})

test_that("RemoveCellGraphs works for Seurat objects", {
  expect_no_error({seur_cleaned <- RemoveCellGraphs(seur_obj)})
  classes <- sapply(seur_cleaned@assays[["mpxCells"]]@cellgraphs, class)
  expect_true(all(classes %in% "NULL"))
})
