library(pixelatorR)
library(dplyr)
library(SeuratObject)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE) %>%
  LoadCellGraphs(cells = colnames(.)[1:2], verbose = FALSE)

cg_list <- CellGraphs(se)[1:2]
set.seed(123)
cg_small_list <- lapply(cg_list, function(cg) {
  nodes <- cg@cellgraph %>%
    pull(name) %>%
    sample(min(5000, length(cg@cellgraph)))
  subset(cg, nodes = nodes)
})
cgl <- CreateCellGraphList(cg_small_list)
cg_small <- cg_small_list[[1]]

test_that("CreateCellGraphList works as expected", {
  expect_s4_class(cgl, "CellGraphList")
  expect_equal(length(cgl), 2)
  expect_equal(names(cgl), names(cg_small_list))
  expect_s4_class(cgl[[1]], "CellGraph")
  expect_equal(length(cgl[1]), 1)
  expect_type(as.list(cgl), "list")
  expect_s4_class(as.list(cgl)[[1]], "CellGraph")
})

test_that("CreateCellGraphList fails when invalid input is provided", {
  expect_error(CreateCellGraphList("Invalid"))
  expect_error(CreateCellGraphList(list(cg_small)))
  expect_error(CreateCellGraphList(list(a = cg_small, a = cg_small)))
  expect_error(CreateCellGraphList(list(a = cg_small, b = "Invalid")))
})

test_that("ComputeLPS.CellGraph stores a matrix as a layer", {
  expect_no_error(cg_lps <- ComputeLPS(cg_small, markers = c("B2M", "CD45")))
  expect_true("lps" %in% Layers(cg_lps))
  lps <- LayerData(cg_lps, layer = "lps")
  expect_equal(nrow(lps), length(Cells(cg_small)))
  expect_equal(ncol(lps), 2)
  expect_equal(colnames(lps), c("B2M", "CD45"))
  expect_equal(rownames(lps), Cells(cg_small))
})

test_that("ComputeLPS.CellGraph stores a vector in meta.data", {
  expect_no_error(cg_lps <- ComputeLPS(cg_small, markers = "B2M", mode = "all", name = "lps_b2m"))
  meta <- CellGraphData(cg_lps, slot = "meta.data")
  expect_true("lps_b2m" %in% colnames(meta))
  expect_equal(nrow(meta), length(Cells(cg_small)))
  expect_false("lps" %in% setdiff(Layers(cg_lps), "counts"))
})

test_that("ComputeLPS.CellGraphList works as expected", {
  expect_no_error(cgl_lps <- ComputeLPS(cgl, markers = "B2M", verbose = FALSE))
  expect_true("lps" %in% Layers(cgl_lps[[1]]))
  expect_true("lps" %in% Layers(cgl_lps[[2]]))
})

test_that("ComputeLPS.CellGraphList intersects missing markers", {
  expect_no_error(
    cgl_lps <- ComputeLPS(cgl, markers = c("B2M", "NotAMarker"), verbose = FALSE)
  )
  expect_equal(colnames(LayerData(cgl_lps[[1]], layer = "lps")), "B2M")
  expect_equal(colnames(LayerData(cgl_lps[[2]], layer = "lps")), "B2M")
})

test_that("ComputeLPS.CellGraphList warns when all markers are missing", {
  expect_warning(
    cgl_skip <- ComputeLPS(cgl[1], markers = "NotAMarker", verbose = FALSE)
  )
  expect_false("lps" %in% setdiff(Layers(cgl_skip[[1]]), "counts"))
})

test_that("ComputeLPS.PNAAssay and Seurat only process loaded cellgraphs", {
  se_one <- se
  cgs <- lapply(CellGraphs(se_one), function(x) NULL)
  cgs[[1]] <- cg_small
  CellGraphs(se_one) <- cgs

  expect_no_error(pna_lps <- ComputeLPS(se_one[["PNA"]], markers = "B2M", verbose = FALSE))
  expect_true("lps" %in% Layers(CellGraphs(pna_lps)[[1]]))
  expect_null(CellGraphs(pna_lps)[[2]])

  expect_no_error(se_lps <- ComputeLPS(se_one, markers = "B2M", verbose = FALSE))
  expect_true("lps" %in% Layers(CellGraphs(se_lps)[[1]]))
  expect_null(CellGraphs(se_lps)[[2]])
})

test_that("ComputeLPS fails when invalid input is provided", {
  expect_error(ComputeLPS("Invalid"))
  expect_error(ComputeLPS(cg_small, markers = "Invalid"))
  expect_error(ComputeLPS(cg_small, mode = "Invalid"))
  expect_error(ComputeLPS(cg_small, method = "Invalid"))
})
