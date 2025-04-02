seur_obj_mpx <- ReadMPX_Seurat(minimal_mpx_pxl_file())
seur_obj_pna <- ReadPNA_Seurat(minimal_pna_pxl_file())

test_that("FSMap works as expected", {
  # Seurat object
  expect_no_error(fs_map <- FSMap(seur_obj_mpx))
  expect_equal(dim(fs_map), c(1, 3))
  expect_equal(names(fs_map), c("id_map", "sample", "pxl_file"))
  expect_no_error(FSMap(seur_obj_mpx) <- FSMap(seur_obj_mpx))

  # CellGraphAssay
  expect_no_error(fs_map_mpx <- FSMap(seur_obj_mpx[["mpxCells"]]))
  expect_equal(dim(fs_map_mpx), c(1, 3))
  expect_equal(names(fs_map_mpx), c("id_map", "sample", "pxl_file"))
  expect_no_error(FSMap(seur_obj_mpx[["mpxCells"]]) <- FSMap(seur_obj_mpx[["mpxCells"]]))

  # PNAAssay
  expect_no_error(fs_map_pna <- FSMap(seur_obj_pna[["PNA"]]))
  expect_equal(dim(fs_map_pna), c(1, 3))
  expect_equal(names(fs_map_pna), c("id_map", "sample", "pxl_file"))
  expect_no_error(FSMap(seur_obj_pna[["PNA"]]) <- FSMap(seur_obj_pna[["PNA"]]))

  # Check that FSMap can be updated
  tmp_path <- fs::file_temp(ext = "pxl")
  fs::file_copy(fs_map$pxl_file, tmp_path)
  expect_no_error({
    FSMap(seur_obj_mpx) <- FSMap(seur_obj_mpx) %>%
      mutate(pxl_file = tmp_path %>% as.character())
  })
  # Check that the update fails if the file is missing
  fs::file_delete(tmp_path)
  expect_error({
    FSMap(seur_obj) <- FSMap(seur_obj) %>%
      mutate(pxl_file = tmp_path %>% as.character())
  })
})

test_that("FSMap fails when invalid input is provided", {
  expect_error(fs_map <- FSMap("Invalid"))
  expect_error(FSMap(seur_obj) <- "Invalid")
})
