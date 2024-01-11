pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
seur_obj[["mpxCells"]] <- KeepLargestComponent(seur_obj[["mpxCells"]])
set.seed(123)
seur_obj <- ComputeLayout(seur_obj)
cg <- seur_obj[["mpxCells"]]@cellgraphs[[colnames(seur_obj)[1]]]

test_that(".add_coordinates_to_tbl_graph works as expected", {

  # Default
  expect_no_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]]))
  expect_true(all(near(cg@cellgraph %>% pull(x) %>% head(n = 2),
               c(0.43121790, 0.01208028))))
  expect_true(all(near(cg@cellgraph %>% pull(y) %>% head(n = 2),
               c(-0.23547434, 0.10504119))))

  # No scaling
  expect_no_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], scale = FALSE))
  expect_true(all(near(cg@cellgraph %>% pull(x) %>% head(n = 2),
                       c(70.476836, 1.974361), tol = 1e-5)))
  expect_true(all(near(cg@cellgraph %>% pull(y) %>% head(n = 2),
                       c(-38.48515, 17.16759), tol = 1e-5)))

  # No scaling
  expect_no_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], keep_aspect_ratio = FALSE))
  expect_true(all(near(cg@cellgraph %>% pull(x) %>% head(n = 2),
                       c(0.44955493, 0.01259398), tol = 1e-5)))
  expect_true(all(near(cg@cellgraph %>% pull(y) %>% head(n = 2),
                       c(-0.2354743, 0.1050412), tol = 1e-5)))
})

test_that(".add_coordinates_to_tbl_graph fails when invalid input is provided", {
  expect_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph("Invalid"), "'cg' must be a 'CellGraph' object")
  expect_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], scale = "Invalid"), "'scale' must be TRUE or FALSE")
  expect_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], keep_aspect_ratio = "Invalid"), "'keep_aspect_ratio' must be TRUE or FALSE")
})
