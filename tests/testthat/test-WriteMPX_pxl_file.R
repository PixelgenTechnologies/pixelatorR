for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  pxl_file <- system.file("extdata/five_cells",
    "five_cells.pxl",
    package = "pixelatorR"
  )
  seur_obj <- ReadMPX_Seurat(pxl_file) %>% LoadCellGraphs()

  seur_obj_merged <- merge(seur_obj, seur_obj, add.cell.ids = c("A", "B")) %>% LoadCellGraphs()

  test_that("WriteMPX_pxl_file works as expected", {
    ## Single sample
    temp_pxl_file <- fs::file_temp(ext = ".pxl")
    expect_no_error({
      WriteMPX_pxl_file(seur_obj, temp_pxl_file)
    })
    expect_true(fs::file_exists(temp_pxl_file))

    # Check content for single sample pxl file
    pxl_files <- unzip(temp_pxl_file, list = TRUE)$Name
    expect_true(all(c(
      "adata.h5ad", "colocalization.parquet", "edgelist.parquet",
      "metadata.json", "polarization.parquet"
    ) %in% pxl_files))

    # Try reloading data
    seur_obj_reloaded <- ReadMPX_Seurat(temp_pxl_file) %>% LoadCellGraphs()
    expect_equal(colnames(seur_obj_reloaded), colnames(seur_obj))
    expect_equal(rownames(seur_obj_reloaded), rownames(seur_obj))
    expect_equal(PolarizationScores(seur_obj_reloaded), PolarizationScores(seur_obj))
    expect_equal(ColocalizationScores(seur_obj_reloaded), ColocalizationScores(seur_obj))

    # Check that the cellgraphs are identical
    cg_list <- CellGraphs(seur_obj)
    cg_list <- lapply(cg_list, function(x) {
      list(cellgraph = x@cellgraph[], counts = x@counts, g_attr = attributes(x@cellgraph))
    })
    cg_list_reloaded <- CellGraphs(seur_obj_reloaded)
    cg_list_reloaded <- lapply(cg_list_reloaded, function(x) {
      list(cellgraph = x@cellgraph[], counts = x@counts, g_attr = attributes(x@cellgraph))
    })
    expect_identical(cg_list, cg_list_reloaded)
    expect_equal(names(cg_list), names(cg_list_reloaded))

    ## Merged sample
    temp_pxl_file <- fs::file_temp(ext = ".pxl")
    expect_no_error({
      WriteMPX_pxl_file(seur_obj_merged, temp_pxl_file)
    })
    expect_true(fs::file_exists(temp_pxl_file))

    # Check content for single sample pxl file
    pxl_files <- unzip(temp_pxl_file, list = TRUE)$Name
    expect_true(all(c(
      "adata.h5ad", "colocalization.parquet", "edgelist.parquet",
      "metadata.json", "polarization.parquet"
    ) %in% pxl_files))

    # Try reloading data
    seur_obj_reloaded <- ReadMPX_Seurat(temp_pxl_file) %>% LoadCellGraphs()
    expect_equal(colnames(seur_obj_reloaded), colnames(seur_obj_merged))
    expect_equal(rownames(seur_obj_reloaded), rownames(seur_obj_merged))
    expect_equal(PolarizationScores(seur_obj_reloaded), PolarizationScores(seur_obj_merged))
    expect_equal(ColocalizationScores(seur_obj_reloaded), ColocalizationScores(seur_obj_merged))

    # Check that the cellgraphs are identical
    cg_list <- CellGraphs(seur_obj_merged)
    cg_list <- lapply(cg_list, function(x) {
      list(cellgraph = x@cellgraph[], counts = x@counts, g_attr = attributes(x@cellgraph)[1:3])
    })
    cg_list_reloaded <- CellGraphs(seur_obj_reloaded)
    cg_list_reloaded <- lapply(cg_list_reloaded, function(x) {
      list(cellgraph = x@cellgraph[], counts = x@counts, g_attr = attributes(x@cellgraph)[1:3])
    })
    expect_identical(cg_list, cg_list_reloaded)
    expect_equal(names(cg_list), names(cg_list_reloaded))
  })
}

if (TRUE) skip("Skipping anndata tests")

## Test anndata
pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
seur_obj <- ReadMPX_Seurat(pxl_file) %>%
  LoadCellGraphs() %>%
  ComputeLayout(dim = 3)
temp_pxl_file <- fs::file_temp(ext = ".pxl")
WriteMPX_pxl_file(seur_obj, temp_pxl_file)

# Unpack adata.h5ad
unzip(temp_pxl_file, files = "adata.h5ad", exdir = fs::path_temp())

# Import anndata python library
anndata <- reticulate::import("anndata")

test_that("anndata file can be loaded with the anndata python library", {
  adata <- anndata$read_h5ad(file.path(fs::path_temp(), "adata.h5ad"))
  expect_equal(rownames(adata$obs), colnames(seur_obj))
  X <- t(adata$X)
  rownames(X) <- rownames(adata$var)
  colnames(X) <- rownames(adata$obs)
  expect_equal(X, LayerData(seur_obj, layer = "counts") %>% as.matrix())
})
