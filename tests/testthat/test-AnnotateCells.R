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
  expect_no_error(seur_anno <- AnnotateCells(
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
    seur_anno$cell_type %>% unname(),
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

  # Test that normalization is preserved
  reference <- read_pbmc_reference()
  seur_sub <-
    seur %>%
    Seurat::NormalizeData(normalization.method = "CLR", margin = 2) %>%
    # Z score scaling
    Seurat::ScaleData() %>%
    Seurat::RunPCA(npcs = 3, features = rownames(seur)) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:3) %>%
    Seurat::FindClusters(resolution = 0.5) %>%
    subset(features = head(rownames(seur), 50))


  # NMF default
  expect_no_error(seur_anno <- AnnotateCells(
    seur_sub,
    reference,
    reference_groups = c("celltype.l1", "celltype.l2"),
    # Summarise annotations by cluster, adding an extra column with majority vote
    summarize_by_column = "seurat_clusters",
    query_assay = "PNA",
    reference_assay = "PNA",
    method = "nmf"
  ))

  expect_equal(
    LayerData(seur_anno, "counts"),
    LayerData(seur_sub, "counts")
  )
  expect_equal(
    LayerData(seur_anno, "data"),
    LayerData(seur_sub, "data")
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

test_that("Reading reference data works as expected", {
  expect_no_error(ref_data <- read_pbmc_reference())

  expect_equal(
    LayerData(ref_data, "counts")[1:10, 1:10],
    new("dgCMatrix",
      i = c(
        0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
        0L, 1L, 2L, 3L, 4L, 5L, 6L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L,
        7L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 0L, 1L, 2L,
        4L, 5L, 6L, 7L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
        0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L,
        7L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 0L, 1L, 2L,
        3L, 4L, 5L, 6L, 7L, 8L, 9L
      ),
      p = c(
        0L, 10L, 19L, 29L, 39L, 48L,
        58L, 67L, 77L, 87L, 97L
      ),
      Dim = c(10L, 10L),
      Dimnames = list(
        c(
          "HLA-ABC", "B2M", "CD11b", "CD11c", "CD18", "CD82", "CD8",
          "TCRab", "HLA-DR", "CD45"
        ),
        c(
          "Resting_PBMC_82c6f1502ceb578c",
          "Resting_PBMC_9524c5af1f13409e", "Resting_PBMC_1a72fb4c1ebd3edf",
          "Resting_PBMC_791355c983aed829", "Resting_PBMC_81415c5542195eff",
          "Resting_PBMC_5ba613cb25de4e5c", "Resting_PBMC_5e08b1ef86b8a4c9",
          "Resting_PBMC_ef963c459c4a8cd6", "Resting_PBMC_89d045aeb08b3401",
          "Resting_PBMC_ce389cbfa22fcb90"
        )
      ),
      x = c(
        1294, 2002, 5, 3,
        132, 1, 11, 7, 2, 1553, 2191, 4690, 3, 6, 295, 5, 4, 6, 1931,
        1368, 3990, 1, 2, 310, 1717, 6, 72, 5, 2098, 17911, 21774, 3,
        12, 18, 1661, 29, 7, 12, 124, 1368, 2357, 2, 146, 25, 1, 10,
        2, 2464, 1482, 2341, 10, 15, 180, 35, 23, 193, 11, 2306, 1252,
        2858, 1, 2, 108, 27, 3390, 97, 1377, 942, 1767, 6, 15, 95, 15,
        477, 6, 1, 976, 2046, 4187, 17, 20, 376, 271, 2872, 95, 25, 2447,
        9771, 11606, 5, 5, 3, 138, 2, 2, 7, 66
      ),
      factors = list()
    )
  )

  expect_equal(
    ref_data$celltype.l1[1:10],
    c(
      Resting_PBMC_82c6f1502ceb578c = "NK", Resting_PBMC_9524c5af1f13409e = "NK",
      Resting_PBMC_1a72fb4c1ebd3edf = "CD4 T", Resting_PBMC_791355c983aed829 = "Platelets",
      Resting_PBMC_81415c5542195eff = "CD4 T", Resting_PBMC_5ba613cb25de4e5c = "CD4 T",
      Resting_PBMC_5e08b1ef86b8a4c9 = "CD8 T", Resting_PBMC_ef963c459c4a8cd6 = "NK",
      Resting_PBMC_89d045aeb08b3401 = "CD8 T", Resting_PBMC_ce389cbfa22fcb90 = "Platelets"
    )
  )
})
