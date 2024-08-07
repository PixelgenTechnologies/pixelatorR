for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  pxl_file <- system.file("extdata/five_cells",
    "five_cells.pxl",
    package = "pixelatorR"
  )
  seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
  seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])

  test_that("RemoveCellGraphs works for MPXAssay objects", {
    expect_no_error({
      cg_assay <- RemoveCellGraphs(seur_obj[["mpxCells"]])
    })
    classes <- sapply(cg_assay@cellgraphs, class)
    expect_true(all(classes %in% "NULL"))
  })

  test_that("RemoveCellGraphs works for Seurat objects", {
    expect_no_error({
      seur_cleaned <- RemoveCellGraphs(seur_obj)
    })
    classes <- sapply(seur_cleaned@assays[["mpxCells"]]@cellgraphs, class)
    expect_true(all(classes %in% "NULL"))
  })
}
