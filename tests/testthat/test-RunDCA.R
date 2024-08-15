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
    targets = "Sample1", reference = "Sample2"
  )))

  expected_result <-
    structure(
      list(
        estimate = c(
          `difference in location` = 0,
          `difference in location` = -5.67168961880486e-05
        ),
        data_type = c("pearson_z", "pearson_z"),
        target = c(
          "Sample1",
          "Sample1"
        ),
        reference = c("Sample2", "Sample2"),
        n1 = c(4L, 4L),
        n2 = c(4L, 4L),
        statistic = c(W = 8, W = 8),
        p = c(1, 1),
        p_adj = c(1, 1),
        conf.low = c(-7.4540803815722, -3.64101112325344),
        conf.high = c(7.4540803815722, 3.6409494101764),
        method = c(
          "Wilcoxon",
          "Wilcoxon"
        ),
        alternative = c("two.sided", "two.sided"),
        marker_1 = c(
          "ACTB",
          "ACTB"
        ),
        marker_2 = c("B2M", "CD102")
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )

  expect_equal(dca_markers[1:2, ], expected_result)

  # cl = 1 should witch to sequential processing
  expect_no_error(suppressWarnings(dca_markers <- RunDCA(colocalization_table_merged,
                                                         contrast_column = "sample",
                                                         targets = "Sample1", reference = "Sample2",
                                                         cl = 1
  )))

  # Colocalization heatmap
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, colors = c("red", "blue")))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, value_col = "statistic"))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_heatmap <- ColocalizationHeatmap(dca_markers, symmetrise = FALSE))
  expect_s3_class(p_heatmap, "pheatmap")
  expect_no_error(p_data <- ColocalizationHeatmap(dca_markers, return_plot_data = TRUE))
  expect_s3_class(p_data, "data.frame")
  expect_equal(dim(p_data), c(80, 80))
})

test_that("RunDCA works as expected on a Seurat object", {
  expect_no_error(suppressWarnings(dca_markers <- RunDCA(seur_merged,
    contrast_column = "sample",
    targets = "Sample1", reference = "Sample2"
  )))

  expected_result <- structure(
    list(
      estimate = c(`difference in location` = 0),
      data_type = "pearson_z",
      target = "Sample1",
      reference = "Sample2",
      n1 = 4L,
      n2 = 4L,
      statistic = c(W = 8),
      p = 1,
      p_adj = 1,
      conf.low = -8.27426094614066,
      conf.high = 8.27426094614066,
      method = "Wilcoxon",
      alternative = "two.sided",
      marker_1 = "ACTB",
      marker_2 = "HLA-ABC"
    ),
    row.names = c(NA, -1L),
    class = c("tbl_df", "tbl", "data.frame")
  )

  expect_equal(dca_markers, expected_result)

  # Automatic selection of targets
  expect_no_error(suppressWarnings(dca_markers <- RunDCA(seur_merged, contrast_column = "sample", reference = "Sample2")))
})


test_that("RunDCA fails with invalid input", {
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged),
    'argument "contrast_column" is missing, with no default'
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample"),
    'argument "reference" is missing, with no default'
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample1"),
    'argument "reference" is missing, with no default'
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "Invalid", targets = "Sample1", reference = "Sample2")
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Invalid", reference = "Sample2")
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample1", reference = "Invalid")
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = c("Sample1", "Sample2"), reference = "Sample1"),
    "targets is invalid\nall targets must be different from reference = 'Sample1'"
  )
  expect_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", group_vars = "sample"),
    "contrast_column = 'sample' cannot be one of group_vars"
  )

  expect_error(
    dpa_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", group_vars = "Invalid"),
    "group_vars is invalid\ngroup_vars must be a character vector with valid column names"
  )
})


if (TRUE) skip("Skipping parallel processing tests")

test_that("RunDCA can be parallelized", {
  # Sequential processing for reference
  expect_no_error(
    dca_markers_seq <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", cl = 2)
  )

  # Using 2 threads. This will be ignored on Windows
  expect_no_error(
    dca_markers_par <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", cl = 2)
  )

  expect_equal(dca_markers_seq, dca_markers_par)

  # Using a cluster object
  cl <- parallel::makeCluster(2)
  expect_no_error(
    dca_markers <- RunDCA(colocalization_table_merged, contrast_column = "sample", targets = "Sample2", reference = "Sample1", cl = cl)
  )
  parallel::stopCluster(cl)
})
