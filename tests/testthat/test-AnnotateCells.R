options(Seurat.object.assay.version = "v3")
seur <- ReadPNA_Seurat(minimal_pna_pxl_file(), overwrite = TRUE, load_proximity_scores = FALSE, verbose = FALSE)
seur <- suppressWarnings(merge(seur, list(seur, seur, seur))) %>%
  Seurat::NormalizeData(normalization.method = "CLR", margin = 2, verbose = FALSE)

set.seed(123)
markers <- rownames(seur)
counts <- Matrix::rowSums(LayerData(seur, "counts"))
sim_ref_counts <- do.call(cbind, rep(list(counts), 30))
colnames(sim_ref_counts) <- paste0("cell", 1:30)
sim_ref_counts["CD4", ] <- c(rep(1e4, 10), rep(0, 20))
sim_ref_counts["CD123", ] <- c(rep(0, 10), rep(1e4, 10), rep(0, 10))
sim_ref_counts["CD16", ] <- c(rep(0, 20), rep(1e4, 10))

reference <- SeuratObject::CreateSeuratObject(
  counts = as(sim_ref_counts, "dgCMatrix"),
  assay = "anot",
  meta.data = data.frame(
    cell_type = c(rep("CD4T", 10), rep("pDC", 10), rep("Mono CD16+", 10)),
    row.names = paste0("cell", 1:30)
  )
) %>%
  Seurat::NormalizeData(normalization.method = "CLR", margin = 2, verbose = FALSE)

test_that("annotate_cells works as expected", {
  # NMF default
  expect_no_error(seur <- AnnotateCells(
    seur,
    reference,
    query_assay = "PNA",
    reference_assay = "anot",
    reference_groups = "cell_type",
    method = "nmf",
    skip_normalization = TRUE,
    verbose = FALSE
  ))
  expect_equal(
    seur$cell_type %>% unname(),
    c(
      "Mono CD16+", "pDC", "CD4T", "CD4T", "CD4T", "Mono CD16+",
      "pDC", "CD4T", "CD4T", "CD4T", "Mono CD16+", "pDC", "CD4T", "CD4T",
      "CD4T", "Mono CD16+", "pDC", "CD4T", "CD4T", "CD4T"
    )
  )

  # NMF threshold
  expect_no_error(seur <- AnnotateCells(
    seur,
    reference,
    query_assay = "PNA",
    reference_assay = "anot",
    reference_groups = "cell_type",
    method = "nmf",
    skip_normalization = TRUE,
    min_prediction_score = 0.5,
    verbose = FALSE
  ))
  expect_equal(
    seur$cell_type %>% unname(),
    c(
      "Mono CD16+", "Unknown", "CD4T", "CD4T", "CD4T", "Mono CD16+",
      "Unknown", "CD4T", "CD4T", "CD4T", "Mono CD16+", "Unknown", "CD4T",
      "CD4T", "CD4T", "Mono CD16+", "Unknown", "CD4T", "CD4T", "CD4T"
    )
  )
})


test_that("annotate_cells fails with invalid input", {
  expect_error(seur <- AnnotateCells("Invalid"))
  expect_error(seur <- AnnotateCells(seur, "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, summarize_by_column = "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, reference_assay = "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, quary_assay = "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, reference_assay = "anot", reference_groups = "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, reference_assay = "anot", reference_groups = "cell_type", method = "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, reference_assay = "anot", reference_groups = "cell_type", method = "nmf", min_prediction_score = "Invalid"))
  expect_error(seur <- AnnotateCells(seur, reference, reference_assay = "anot", reference_groups = "cell_type", method = "nmf", skip_normalization = "Invalid"))
})
