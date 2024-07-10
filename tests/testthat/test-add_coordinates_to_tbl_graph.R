options(Seurat.object.assay.version = "v3")

pxl_file <- system.file("extdata/five_cells",
                        "five_cells.pxl",
                        package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])

set.seed(123)
seur_obj <- ComputeLayout(seur_obj, layout_method = "pmds")
cg <- CellGraphs(seur_obj)[[1]]

test_that(".add_coordinates_to_tbl_graph works as expected", {

  # Default
  expect_no_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]]))
  expect_equal(
    cg@cellgraph %>% select(x, y) %>% as_tibble() %>% head(n = 10),
    structure(
      list(
        x = c(
          -0.319142999539692,
          -0.258939989825134,
          -0.142907851240755,
          0.329872509477288,-0.260545413452768,
          -0.302551946973838,
          0.65282972071568,
          -0.0243226680565593,
          0.665851086396988,
          0.411071497198462
        ),
        y = c(
          -0.0240512472495484,
          -0.510750890649215,
          -0.0275715484388217,
          0.0100691569069965,
          0.242934495621249,
          -0.168168698028061,
          -0.0640596552246973,
          -0.00822271995019247,-0.255294040606221,
          -0.205507780956781
        )
      ),
      row.names = c(NA,-10L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  # No scaling
  expect_no_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], scale = FALSE))
  expect_equal(
    cg@cellgraph %>% select(x, y) %>% as_tibble() %>% head(n = 10),
    structure(
      list(
        x = c(
          -58.5281066499656,
          -47.4873876672378,
          -26.208082178209,
          60.4958073447489,
          -47.7818086804891,
          -55.4854489842481,
          119.723407918783,-4.46057006445151,
          122.111415427715,
          75.387009783996
        ),
        y = c(
          -4.41079380126333,-93.6673611597246,
          -5.05638704235889,
          1.8465975759391,
          44.5521164154744,-30.8407062275168,
          -11.7479949062295,
          -1.50797677182892,
          -46.8187516481721,-37.6883758646892
        )
      ),
      row.names = c(NA,-10L),
      class = c("tbl_df",
                "tbl", "data.frame")
    )
  )

  # No scaling
  expect_no_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], keep_aspect_ratio = FALSE))
  expect_equal(
    cg@cellgraph %>% select(x, y) %>% as_tibble() %>% head(n = 10),
    structure(
      list(
        x = c(
          -0.319142999539692,
          -0.258939989825134,-0.142907851240755,
          0.329872509477288,
          -0.260545413452768,
          -0.302551946973838,
          0.65282972071568,
          -0.0243226680565593,
          0.665851086396988,
          0.411071497198462
        ),
        y = c(
          -0.0297533213167315,
          -0.631839804589686,
          -0.0341082161514718,
          0.0124563544557166,
          0.300529449975732,
          -0.208038163506852,
          -0.07924692992274,-0.0101721638897612,
          -0.315819198131033,
          -0.254229602999503
        )
      ),
      row.names = c(NA,-10L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that(".add_coordinates_to_tbl_graph fails when invalid input is provided", {
  expect_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph("Invalid"), "'cg' must be a 'CellGraph' object")
  expect_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], scale = "Invalid"), "'scale' must be TRUE or FALSE")
  expect_error(cg <- pixelatorR:::.add_coordinates_to_tbl_graph(cg, layout_coordinates = cg@layout[["pmds"]], keep_aspect_ratio = "Invalid"), "'keep_aspect_ratio' must be TRUE or FALSE")
})
