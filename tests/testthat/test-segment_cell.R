options(Seurat.object.assay.version = "v5")
se <- ReadPNA_Seurat(minimal_pna_pxl_file()) %>%
  LoadCellGraphs(cells = colnames(.)[2], verbose = FALSE)
cg <- CellGraphs(se)[[2]]

# Get cell type weights
se$cell_type <- c("Mono", "pDC", "CD4T", "CD4T", "CD4T")
set.seed(7331)
w <- cc_protein_weights(
  se,
  group_by = "cell_type",
  population_1 = "Mono", population_2 = "CD4T",
  mode = "cell_abundance",
  show_plot = FALSE
)

test_that("segment_cell returns expected compartments", {
  expect_no_error(
    out <- segment_cell(
      cg,
      w = w,
      k = 2L,
      verbose = FALSE
    )
  )

  comp <- out@cellgraph %>% pull(compartment)
  expect_true(all(comp %in% c("Mono", "CD4T", "interface", "other")))
  expect_true(any(comp %in% c("Mono", "CD4T")))
  expect_false(any(grepl("^intra_", comp)))
})

test_that("segment_cell without interface detection has no interface labels", {
  out <- segment_cell(
    cg,
    w = w,
    k = 2L,
    detect_interface = FALSE,
    verbose = FALSE
  )

  comp <- out@cellgraph %>% pull(compartment)
  expect_false(any(comp == "interface"))
  expect_true(all(comp %in% c("Mono", "CD4T", "other")))
})

test_that("segment_cell handles empty interface without error", {
  # Very high expansion tends to absorb all segmented nodes into the interface,
  # leaving no crossing edges for residual interface detection in downstream logic.
  expect_no_error(
    out <- segment_cell(
      cg,
      w = w,
      k = 2L,
      k_interface_expansion = 4L,
      verbose = FALSE
    )
  )

  comp <- out@cellgraph %>% pull(compartment)
  expect_true(all(comp %in% c("Mono", "CD4T", "interface", "other")))
})

test_that("segment_cell applies component filtering to final annotation", {
  out_lcc <- segment_cell(
    cg,
    w = w,
    k = 2L,
    detect_interface = FALSE,
    keep_largest_comp = TRUE,
    verbose = FALSE
  )

  out_min_size <- segment_cell(
    cg,
    w = w,
    k = 2L,
    detect_interface = FALSE,
    keep_largest_comp = FALSE,
    min_comp_size = 1L,
    verbose = FALSE
  )

  comp_lcc <- out_lcc@cellgraph %>% pull(compartment)
  comp_min_size <- out_min_size@cellgraph %>% pull(compartment)

  n_non_other_lcc <- sum(comp_lcc %in% c("Mono", "CD4T"))
  n_non_other_min_size <- sum(comp_min_size %in% c("Mono", "CD4T"))

  expect_lte(n_non_other_lcc, n_non_other_min_size)
  expect_gt(sum(comp_lcc == "other"), sum(comp_min_size == "other"))
})

test_that("segment_cell fails fast when w is NULL", {
  expect_error(
    segment_cell(cg, w = NULL),
    regexp = "w.*must be provided"
  )
})

test_that("segment_cell fails with invalid input", {
  expect_error(segment_cell("Invalid", w = w))
  expect_error(segment_cell(cg, w = "Invalid"))
  expect_error(segment_cell(cg, w = w, k = "Invalid"))
  expect_error(segment_cell(cg, w = w, detect_interface = "Invalid"))
  expect_error(segment_cell(cg, w = w, k_interface_expansion = "Invalid"))
  expect_error(segment_cell(cg, w = w, keep_largest_comp = "Invalid"))
  expect_error(segment_cell(cg, w = w, min_comp_size = "Invalid"))
  expect_error(segment_cell(cg, w = w, verbose = "Invalid"))
})
