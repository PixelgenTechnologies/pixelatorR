options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/mock_data", "mock_mpx_data.pxl", package = "pixelatorR")
el <- ReadMPX_arrow_edgelist(pxl_file, overwrite = TRUE)

test_that("edgelist_to_simple_Anode_graph works as expected", {

  # data.frame
  el_tbl_df <- as_tibble(el)
  expect_no_error(anode <- edgelist_to_simple_Anode_graph(el_tbl_df, components = "RCVCMP0000217"))

  # FileSystemDataset
  expect_no_error(anode <- edgelist_to_simple_Anode_graph(el, components = "RCVCMP0000217"))
})

test_that("edgelist_to_simple_Anode_graph fails when invalid input is provided", {

  # data.frame
  el_tbl_df <- as_tibble(el)
  expect_error(anode <- edgelist_to_simple_Anode_graph(el_tbl_df, components = "Invalid"))

  # FileSystemDataset
  expect_error(anode <- edgelist_to_simple_Anode_graph(el, components = "Invalid"))

})

test_that("edgelist_to_simple_bipart_graph works as expected", {

  # data.frame
  el_tbl_df <- as_tibble(el) %>% filter(component == "RCVCMP0000217")
  expect_no_error(anode <- edgelist_to_simple_bipart_graph(el_tbl_df))

})
