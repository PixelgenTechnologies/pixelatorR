pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  test_that("ReadPNA_Seurat works as expected", {
    # default
    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, verbose = FALSE)))
    expect_s4_class(seur_obj, "Seurat")
    expect_equal(dim(seur_obj), c(158, 5))
    expect_s4_class(seur_obj[["PNA"]], ifelse(assay_version == "v3", "PNAAssay", "PNAAssay5"))

    # no proximity scores
    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, load_proximity_scores = FALSE, verbose = FALSE)))
    expect_true(length(ProximityScores(seur_obj)) == 0)

    # new assay name
    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, assay = "myassay", verbose = FALSE)))
    expect_true("myassay" %in% SeuratObject::Assays(seur_obj))

    expect_no_error(suppressWarnings(seur_obj <- ReadPNA_Seurat(pxl_file, return_pna_assay = FALSE, verbose = FALSE)))
    expect_s4_class(seur_obj[["PNA"]], ifelse(assay_version == "v3", "Assay", "Assay5"))
  })

  test_that("ReadPNA_Seurat fails with invalid input", {
    expect_error(ReadPNA_Seurat("Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, verbose = "Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, return_pna_assay = "Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, load_proximity_scores = "Invalid"))
    expect_error(ReadPNA_Seurat(pxl_file, assay = FALSE))
  })
}
