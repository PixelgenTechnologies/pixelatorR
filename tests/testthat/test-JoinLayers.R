options(Seurat.object.assay.version = "v5")

pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file)

# Merge Seurat objects
seur_obj_merged <- merge(seur_obj, seur_obj, add.cell.ids = c("A", "B"))
cg_assay5 <- seur_obj_merged[["mpxCells"]]

test_that("JoinLayers works as expected", {

  # Seurat object
  expect_no_error(seur_obj_merged <- JoinLayers(seur_obj_merged))
  expect_true(is(seur_obj_merged[["mpxCells"]], "CellGraphAssay5"))
  expect_equal(seur_obj_merged[["mpxCells"]]@layers %>% length(), 1)

  # CellGraphAssay5 object
  expect_no_error(cg_assay5 <- JoinLayers(cg_assay5))
  expect_true(is(cg_assay5, "CellGraphAssay5"))
  expect_equal(cg_assay5@layers %>% length(), 1)
})
