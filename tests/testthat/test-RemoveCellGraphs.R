mpx_pxl_file <- minimal_mpx_pxl_file()
pna_pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  seur_obj_mpx <- ReadMPX_Seurat(mpx_pxl_file)
  seur_obj_mpx <- LoadCellGraphs(seur_obj_mpx, cells = colnames(seur_obj_mpx)[1], verbose = FALSE)

  seur_obj_pna <- ReadPNA_Seurat(pna_pxl_file)
  seur_obj_pna <- LoadCellGraphs(seur_obj_pna, cells = colnames(seur_obj_pna)[1], verbose = FALSE)

  test_that("RemoveCellGraphs works for MPXAssay objects", {
    expect_no_error({
      cg_assay <- RemoveCellGraphs(seur_obj_mpx[["mpxCells"]])
    })
    expect_no_error({
      pna_assay <- RemoveCellGraphs(seur_obj_pna[["PNA"]])
    })
    classes_mpx <- sapply(cg_assay@cellgraphs, class)
    classes_pna <- sapply(pna_assay@cellgraphs, class)
    expect_true(all(classes_mpx %in% "NULL"))
    expect_true(all(classes_pna %in% "NULL"))
  })

  test_that("RemoveCellGraphs works for Seurat objects", {
    expect_no_error({
      seur_mpx_cleaned <- RemoveCellGraphs(seur_obj_mpx)
    })
    expect_no_error({
      seur_pna_cleaned <- RemoveCellGraphs(seur_obj_pna)
    })
    classes_mpx <- sapply(seur_mpx_cleaned@assays[["mpxCells"]]@cellgraphs, class)
    expect_true(all(classes_mpx %in% "NULL"))
    classes_pna <- sapply(seur_pna_cleaned@assays[["PNA"]]@cellgraphs, class)
    expect_true(all(classes_pna %in% "NULL"))
  })
}
