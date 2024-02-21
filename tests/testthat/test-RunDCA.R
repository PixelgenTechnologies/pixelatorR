library(dplyr)
options(pixelatorR.arrow_outdir = tempdir())

pxl_file <- system.file("extdata/PBMC_10_cells",
                        "Sample01_test.pxl",
                        package = "pixelatorR")

# Load colocalization scores
colocalization_table1 <- colocalization_table2 <- ReadMPX_colocalization(pxl_file)
colocalization_table1$sample <- "Sample1"
colocalization_table2$sample <- "Sample2"
colocalization_table_merged <-  bind_rows(colocalization_table1, colocalization_table2) %>%
  filter(marker_1 %in% c("ACTB", "HLA-ABC"))

# Seurat objects
seur1 <- seur2 <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
seur1$sample <- "Sample1"
seur2$sample <- "Sample2"
seur_merged <- merge(seur1, seur2, add.cell.ids = c("Sample1", "Sample2"))
seur_merged <- subset(seur_merged, features = c("ACTB", "HLA-ABC"))

test_that("RunDCA works as expected on a data.frame and that ColocalizationHeatmap works on the output", {

  expect_no_error(suppressWarnings(dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample",
                                                         target = "Sample1", reference = "Sample2")))

  expect_true(all(dca_markers$data_type == "pearson_z"))
  expect_true(all(dca_markers$target == "Sample1"))
  expect_true(all(dca_markers$reference == "Sample2"))
  expect_true(all(dca_markers$n1[1:4] == c(9, 5, 7, 10)))
  expect_true(all(dca_markers$n2[1:4] == c(9, 5, 7, 10)))
  expect_true(all(dca_markers$statistic[1:4] == c(40.5, 12.5, 24.5 , 50)))
  expect_true(all(dca_markers$method == "Wilcoxon"))
  expect_true(all(dca_markers$marker_1[1:4] == c("ACTB", "ACTB", "ACTB", "ACTB")))
  expect_true(all(dca_markers$marker_2[1:4] == c("CD11c", "CD14", "CD19", "CD20")))

  # Colocalization heatmap
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, colors = c("red", "blue")))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_data <- ColocalizationHeatmap(dca_markers, return_plot_data = TRUE))
  expect_s3_class(p_data, "data.frame")
  expect_equal(dim(p_data), c(26, 26))
})

test_that("RunDCA works as expected on a Seurat object", {

  expect_no_error(suppressWarnings(dca_markers <- RunDCA(seur_merged, contrast_column = "sample",
                                                         target = "Sample1", reference = "Sample2")))

  expect_true(all(dca_markers$data_type == "pearson_z"))
  expect_true(all(dca_markers$target == "Sample1"))
  expect_true(all(dca_markers$reference == "Sample2"))
  expect_true(all(dca_markers$n1 == 10))
  expect_true(all(dca_markers$n2 == 10))
  expect_true(all(dca_markers$statistic == 50))
  expect_true(all(dca_markers$method == "Wilcoxon"))
  expect_true(all(dca_markers$marker_1 == "ACTB"))
  expect_true(all(dca_markers$marker_2 == "HLA-ABC"))
})


test_that("RunDCA fails with invalid input",  {
  expect_error(dca_markers <- RunDCA(colocalization_table_merged),
               'argument "contrast_column" is missing, with no default')
  expect_error(dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample"),
               'argument "target" is missing, with no default')
  expect_error(dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", target = "Sample1"),
               'argument "reference" is missing, with no default')
  expect_error(dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "Invalid", target = "Sample1", reference = "Sample2"),
               "'contrast_column' must be a valid column name")
  expect_error(dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", target = "Invalid", reference = "Sample2"),
               "'target' must be present in 'contrast_column' column")
  expect_error(dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", target = "Sample1", reference = "Invalid"),
               "'reference' must be present in 'contrast_column' column")
})
