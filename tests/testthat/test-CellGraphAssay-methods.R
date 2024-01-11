pxl_file <- system.file("extdata/PBMC_10_cells", "Sample01_test.pxl", package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
cg_assay <- seur_obj[["mpxCells"]]

# CellGraphs method
test_that("CellGraphs.CellGraphAssay getter/setter works as expected", {
  cg_list <- CellGraphs(cg_assay)
  expect_type(cg_list, "list")
  expect_equal(cg_list %>% length(), 10)
  CellGraphs(cg_assay) <- cg_assay@cellgraphs
  cg_list <- cg_assay@cellgraphs
  expect_equal(cg_list %>% length(), 10)
  CellGraphs(cg_assay) <- NULL
})

test_that("CellGraphs.CellGraphAssay getter/setter fails when invalid input is provided", {
  expect_error(CellGraphs("Invalid input"))
  expect_error(CellGraphs(cg_assay) <- "Invalid input", "Invalid class 'character'")
  expect_error(CellGraphs(cg_assay) <- setNames(cg_assay@cellgraphs, nm = paste0(names(cg_assay@cellgraphs), "invalid")))
  cgs <- cg_assay@cellgraphs
  cgs[[1]] <- "Invalid"
  expect_error({CellGraphs(cg_assay) <- cgs})
})

# RenameCells method
test_that("RenameCells.CellGraphAssay method works as expected", {
  cg_assay_renamed <- RenameCells(cg_assay, new.names = paste0("A_", colnames(cg_assay)))
  expect_s4_class(cg_assay_renamed, "CellGraphAssay")
  expect_equal(colnames(cg_assay_renamed), paste0("A_", colnames(cg_assay)))
})

test_that("RenameCells.CellGraphAssay method fails when invalid input is provided", {
  expect_error(RenameCells(cg_assay, new.names = "Invalid"), "'new.names' must be a character vector with the same length as the number of cells present in 'object'")
})

# subset method
test_that("subset.CellGraphAssay works as expected", {
  cg_assay_subset <- subset(cg_assay, cells = colnames(cg_assay)[1:2])
  expect_equal(ncol(cg_assay_subset), 2)
  expect_equal(colnames(cg_assay_subset), c("RCVCMP0000000", "RCVCMP0000002"))
})

# merge method
cg_assay <- seur_obj[["mpxCells"]]
test_that("merge.CellGraphAssay works as expected", {
  cg_assay_merged <- merge(cg_assay, y = cg_assay)
  expect_equal(ncol(cg_assay_merged), 20)
  expect_equal(length(CellGraphs(cg_assay_merged)), 20)
  cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay))
  expect_equal(ncol(cg_assay_merged), 30)
  expect_equal(length(CellGraphs(cg_assay_merged)), 30)
  cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay))
  expect_equal(colnames(cg_assay_merged), c(paste0(colnames(cg_assay), "_1"),
                                            paste0(colnames(cg_assay), "_2"),
                                            paste0(colnames(cg_assay), "_3")))
})

test_that("merge.CellGraphAssay fails when invalid input is provided", {
  expect_error({cg_assay_merged <- merge(cg_assay, y = "Invalid")}, "'y' must be a 'CellGraphAssay' object or a list of 'CellGraphAssay' objects")
  expect_error({cg_assay_merged <- merge(cg_assay, y = list(cg_assay, "Invalid"))}, "Element 2 in 'y' is not a 'CellGraphAssay'")
})

# Show method
test_that("show.CellGraphAssay works as expected", {
  expect_no_error(capture.output(show(cg_assay)))
})

# PolarizationScores getter/setter method
test_that("PolarizationScores.CellGraphAssay works as expected", {

  # Getter
  expect_no_error(pol <- PolarizationScores(cg_assay))
  expect_s3_class(pol, "tbl_df")
  expect_equal(names(pol), c("morans_i", "morans_p_value", "morans_p_adjusted", "morans_z", "marker", "component"))
  expect_equal(dim(pol), c(225, 6))

  # Setter
  expect_no_error(PolarizationScores(cg_assay) <- PolarizationScores(cg_assay))
  expect_no_error(pol <- PolarizationScores(cg_assay))
  expect_s3_class(pol, "tbl_df")
  expect_equal(names(pol), c("morans_i", "morans_p_value", "morans_p_adjusted", "morans_z", "marker", "component"))
  expect_equal(dim(pol), c(225, 6))
})

test_that("PolarizationScores.CellGraphAssay fails when invalid input is provided", {

  # Getter
  expect_error(pol <- PolarizationScores("invalid"), "no applicable method for")

  # Setter
  expect_error(PolarizationScores(cg_assay) <- "invalid", "'polarization' must be a non-empty 'tbl_df' object")
})


# ColocalizationScores getter/setter method
test_that("ColocalizationScores.CellGraphAssay works as expected", {

  # Getter
  expect_no_error(coloc <- ColocalizationScores(cg_assay))
  expect_s3_class(coloc, "tbl_df")
  expect_equal(names(coloc),
               c("marker_1", "marker_2", "pearson",
                 "pearson_mean", "pearson_stdev", "pearson_z",
                 "pearson_p_value", "pearson_p_value_adjusted",
                 "jaccard", "jaccard_mean", "jaccard_stdev",
                 "jaccard_z", "jaccard_p_value", "jaccard_p_value_adjusted",
                 "component"))
  expect_equal(dim(coloc), c(2680, 15))

  # Setter
  expect_no_error(ColocalizationScores(cg_assay) <- ColocalizationScores(cg_assay))
  expect_no_error(coloc <- ColocalizationScores(cg_assay))
  expect_s3_class(coloc, "tbl_df")
  expect_equal(names(coloc),
               c("marker_1", "marker_2", "pearson",
                 "pearson_mean", "pearson_stdev", "pearson_z",
                 "pearson_p_value", "pearson_p_value_adjusted",
                 "jaccard", "jaccard_mean", "jaccard_stdev",
                 "jaccard_z", "jaccard_p_value", "jaccard_p_value_adjusted",
                 "component"))
  expect_equal(dim(coloc), c(2680, 15))
})

test_that("ColocalizationScores.CellGraphAssay fails when invalid input is provided", {

  # Getter
  expect_error(coloc <- ColocalizationScores("invalid"), "no applicable method for")

  # Setter
  expect_error(ColocalizationScores(cg_assay) <- "invalid", "'colocalization' must be a non-empty 'tbl_df' object")
})


# ArrowData getter/setter method
test_that("ArrowData.CellGraphAssay works as expected", {

  # Getter
  expect_no_error(arrow_data <- ArrowData(cg_assay))
  expect_s3_class(arrow_data, "FileSystemDataset")
  expect_equal(names(arrow_data),
               c("upia", "upib", "umi", "marker", "sequence",
                 "count", "umi_unique_count", "upi_unique_count",
                 "component", "sample"))
  expect_equal(dim(arrow_data), c(64750, 10))

  # Setter
  expect_no_error(ArrowData(cg_assay) <- ArrowData(cg_assay))
  expect_no_error(arrow_data <- ArrowData(cg_assay))
  expect_s3_class(arrow_data, "FileSystemDataset")
  expect_equal(names(arrow_data),
               c("upia", "upib", "umi", "marker", "sequence",
                 "count", "umi_unique_count", "upi_unique_count",
                 "component", "sample"))
  expect_equal(dim(arrow_data), c(64750, 10))
})

test_that("ArrowData.CellGraphAssay fails when invalid input is provided", {

  # Getter
  expect_error(arrow_data <- ArrowData("invalid"), "no applicable method for")

  # Setter
  expect_error(ArrowData(cg_assay) <- "invalid", "Invalid class 'character'")
})


# ArrowDir getter/setter method
test_that("ArrowDir.CellGraphAssay works as expected", {

  # Getter
  expect_no_error(arrow_dir <- ArrowDir(cg_assay))
  expect_type(arrow_dir, "character")
  expect_length(arrow_dir, 1)

  # Setter
  expect_no_error(ArrowDir(cg_assay) <- ArrowDir(cg_assay))
  expect_no_error(arrow_dir <- ArrowDir(cg_assay))
  expect_type(arrow_dir, "character")
  expect_length(arrow_dir, 1)
})

test_that("ArrowDir.CellGraphAssay fails when invalid input is provided", {

  # Getter
  expect_error(arrow_dir <- ArrowDir("invalid"), "no applicable method for")

  # Setter missing directory
  expect_error(ArrowDir(cg_assay) <- "invalid", "invalid doesn't exist")

  # Setter invalid inout
  expect_error(ArrowDir(cg_assay) <- 1, "'value' must be a non-empty character vector")
})
