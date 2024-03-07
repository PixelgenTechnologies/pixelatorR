pxl_file <- system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR")
seur <- ReadMPX_Seurat(pxl_file, return_cellgraphassay = T, overwrite = T)
cg_assay <- seur[["mpxCells"]]

test_that("RestoreArrowConnection works as expected", {

  # CellGraphAssay
  expect_no_error(anode <- RestoreArrowConnection(cg_assay))

  # Seurat
  expect_no_error(seur <- RestoreArrowConnection(seur))
})

test_that("edgelist_to_simple_Anode_graph fails when invalid input is provided", {

  # data.frame
  el_tbl_df <- as_tibble(cg_assay@arrow_data)
  expect_error(anode <- edgelist_to_simple_Anode_graph(el_tbl_df, components = "Invalid"))

  # FileSystemDataset
  expect_error(anode <- edgelist_to_simple_Anode_graph(el, components = "Invalid"))

})

test_that("edgelist_to_simple_bipart_graph works as expected", {

  # data.frame
  el_tbl_df <- as_tibble(cg_assay@arrow_data) %>% filter(component == "RCVCMP0000217")
  expect_no_error(anode <- edgelist_to_simple_bipart_graph(el_tbl_df))

})
