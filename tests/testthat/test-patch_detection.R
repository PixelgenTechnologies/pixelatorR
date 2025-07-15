se <- ReadPNA_Seurat(minimal_pna_pxl_file(), load_proximity_scores = FALSE, verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1], add_layouts = TRUE, verbose = FALSE)

test_that("patch_detection works as expected", {
  expect_no_error(cg_patch <-
    patch_detection(
      CellGraphs(se)[[1]],
      patch_markers = c("CD58"),
      verbose = FALSE
    ))
  expect_no_error(cg_patch@cellgraph %>% pull(patch))
  expect_equal(
    table(cg_patch@cellgraph %>% pull(patch)),
    structure(
      c(`1` = 41588L, `2` = 1267L, `3` = 688L),
      dim = 3L,
      dimnames = structure(list(c("0", "1", "2")), names = ""),
      class = "table"
    )
  )
  expect_no_error(
    cg_patch <-
      patch_detection(
        CellGraphs(se)[[1]],
        patch_markers = c("CD58"),
        patch_nodes_threshold = 10,
        verbose = FALSE,
        method = "local_G"
      )
  )
})

test_that("patch_detection fails with invalid input", {
  expect_error(patch_detection("Invalid"))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", host_markers = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", k = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", leiden_resolution = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", patch_nodes_threshold = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", prune_patch_edge = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", leiden_refinement = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", method = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", expand_contract_scale_factor = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", local_g_pval_threshold = "Invalid", verbose = FALSE))
  expect_error(patch_detection(CellGraphs(se)[[1]], patch_markers = "CD58", seed = "Invalid", verbose = FALSE))
})
