options(Seurat.object.assay.version = "v3")

pxl_file <- system.file("extdata/five_cells", "five_cells.pxl", package = "pixelatorR")
seur_obj <- ReadMPX_Seurat(pxl_file, overwrite = TRUE)
cg_assay <- seur_obj[["mpxCells"]]

# CellGraphs method
test_that("CellGraphs.CellGraphAssay getter/setter works as expected", {
  cg_list <- CellGraphs(cg_assay)
  expect_type(cg_list, "list")
  expect_equal(cg_list %>% length(), 5)
  CellGraphs(cg_assay) <- cg_assay@cellgraphs
  cg_list <- cg_assay@cellgraphs
  expect_equal(cg_list %>% length(), 5)
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
  expect_error(RenameCells(cg_assay, new.names = "Invalid"),
               "'new.names' must be a character vector where ")
})

# subset method
test_that("subset.CellGraphAssay works as expected", {
  cg_assay_subset <- subset(cg_assay, cells = colnames(cg_assay)[1:2])
  expect_equal(ncol(cg_assay_subset), 2)
  expect_equal(colnames(cg_assay_subset), c("RCVCMP0000217", "RCVCMP0000118"))
})

# merge method
cg_assay <- seur_obj[["mpxCells"]]
test_that("merge.CellGraphAssay works as expected", {
  expect_no_error(cg_assay_merged <- merge(cg_assay, y = cg_assay, add.cell.ids = c("A", "B")))
  expect_equal(ncol(cg_assay_merged), 10)
  expect_equal(length(CellGraphs(cg_assay_merged)), 10)
  expect_no_error(cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay), add.cell.ids = c("A", "B", "C")))
  expect_equal(ncol(cg_assay_merged), 15)
  expect_equal(length(CellGraphs(cg_assay_merged)), 15)
  expect_no_error(cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay), add.cell.ids = c("A", "B", "C")))
  expect_equal(colnames(cg_assay_merged), c(paste0("A_", colnames(cg_assay)),
                                            paste0("B_", colnames(cg_assay)),
                                            paste0("C_", colnames(cg_assay))))
  expect_no_error({cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay), add.cell.ids = c("A", "B", "C"))})
  expect_no_error({cg_assay_double_merged <- merge(cg_assay_merged, cg_assay_merged, add.cell.ids = c("A", "B"))})
})

test_that("merge.CellGraphAssay fails when invalid input is provided", {
  expect_error({cg_assay_merged <- merge(cg_assay, y = "Invalid")}, "'y' must be a 'CellGraphAssay")
  expect_error({cg_assay_merged <- merge(cg_assay, y = list(cg_assay, "Invalid"))}, "Element 2 in 'y' is not a")
  expect_no_error({cg_assay_merged <- merge(cg_assay, y = list(cg_assay, cg_assay), add.cell.ids = c("A", "B", "C"))})
  expect_no_error({cg_assay_double_merged <- merge(cg_assay_merged, cg_assay_merged, add.cell.ids = c("A", "B"))})
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
  expect_equal(dim(pol), c(380, 6))

  # Setter
  expect_no_error(PolarizationScores(cg_assay) <- PolarizationScores(cg_assay))
  expect_no_error(pol <- PolarizationScores(cg_assay))
  expect_s3_class(pol, "tbl_df")
  expect_equal(names(pol), c("morans_i", "morans_p_value", "morans_p_adjusted", "morans_z", "marker", "component"))
  expect_equal(dim(pol), c(380, 6))
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
  expect_equal(dim(coloc), c(14649, 15))

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
  expect_equal(dim(coloc), c(14649, 15))
})

test_that("ColocalizationScores.CellGraphAssay fails when invalid input is provided", {

  # Getter
  expect_error(coloc <- ColocalizationScores("invalid"), "no applicable method for")

  # Setter
  expect_error(ColocalizationScores(cg_assay) <- "invalid", "'colocalization' must be a non-empty 'tbl_df' object")
})
