library(dplyr)
options(Seurat.object.assay.version = "v3")

pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)

# Load colocalization scores
colocalization_table1 <- colocalization_table2 <- ReadMPX_colocalization(pxl_file)
colocalization_table1$sample <- "Sample1"
colocalization_table2$sample <- "Sample2"
colocalization_table_merged <- bind_rows(colocalization_table1, colocalization_table2) %>%
  filter(marker_1 %in% c("ACTB", "HLA-ABC"))

# Seurat objects
seur1 <- seur2 <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur1$sample <- "Sample1"
seur2$sample <- "Sample2"
seur_merged <- merge(seur1, seur2, add.cell.ids = c("Sample1", "Sample2"))
seur_merged <- subset(seur_merged, features = c("ACTB", "HLA-ABC"))

test_that("RunDCA works as expected on a data.frame and that ColocalizationHeatmap works on the output", {
  expect_no_error(suppressWarnings(dca_markers <- RunDCA(colocalization_table_merged,
    contrast_column = "sample",
    target = "Sample1", reference = "Sample2"
  )))

  expect_true(all(dca_markers$data_type == "pearson_z"))
  expect_true(all(dca_markers$target == "Sample1"))
  expect_true(all(dca_markers$reference == "Sample2"))
  expect_true(all(dca_markers$n1[1:4] == c(4, 4, 4, 3)))
  expect_true(all(dca_markers$n2[1:4] == c(4, 4, 4, 3)))
  expect_true(all(dca_markers$statistic[1:4] == c(8, 8, 8, 4.5)))
  expect_true(all(dca_markers$method == "Wilcoxon"))
  expect_true(all(dca_markers$marker_1[1:4] == c("ACTB", "ACTB", "ACTB", "ACTB")))
  expect_true(all(dca_markers$marker_2[1:4] == c("B2M", "CD102", "CD11a", "CD11b")))

  # Colocalization heatmap
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, colors = c("red", "blue")))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, value_col = "statistic"))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, symmetrise = F))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_data <- ColocalizationHeatmap(dca_markers, return_plot_data = TRUE))
  expect_s3_class(p_data, "data.frame")
  expect_equal(dim(p_data), c(80, 80))
})

test_that("RunDCA works as expected on a Seurat object", {
  expect_no_error(suppressWarnings(dca_markers <- RunDCA(seur_merged,
    contrast_column = "sample",
    target = "Sample1", reference = "Sample2"
  )))

  expect_true(all(dca_markers$data_type == "pearson_z"))
  expect_true(all(dca_markers$target == "Sample1"))
  expect_true(all(dca_markers$reference == "Sample2"))
  expect_true(all(dca_markers$n1 == 4))
  expect_true(all(dca_markers$n2 == 4))
  expect_true(all(dca_markers$statistic == 8))
  expect_true(all(dca_markers$method == "Wilcoxon"))
  expect_true(all(dca_markers$marker_1 == "ACTB"))
  expect_true(all(dca_markers$marker_2 == "HLA-ABC"))
})


test_that("RunDCA fails with invalid input", {
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged),
    'argument "contrast_column" is missing, with no default'
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample"),
    'argument "target" is missing, with no default'
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", target = "Sample1"),
    'argument "reference" is missing, with no default'
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "Invalid", target = "Sample1", reference = "Sample2"),
    "'contrast_column' must be a valid column name"
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", target = "Invalid", reference = "Sample2"),
    "'target' must be present in 'contrast_column' column"
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", target = "Sample1", reference = "Invalid"),
    "'reference' must be present in 'contrast_column' column"
  )
})
