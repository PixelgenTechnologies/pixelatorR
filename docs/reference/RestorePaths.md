# Restore PXL file paths

Updates the PXL file paths in an object with the paths in the specified
directory.

## Usage

``` r
RestorePaths(object, ...)

# S3 method for class 'MPXAssay'
RestorePaths(object, pxl_files_dir, verbose = TRUE, ...)

# S3 method for class 'CellGraphAssay'
RestorePaths(object, pxl_files_dir, verbose = TRUE, ...)

# S3 method for class 'CellGraphAssay5'
RestorePaths(object, pxl_files_dir, verbose = TRUE, ...)

# S3 method for class 'PNAAssay'
RestorePaths(object, pxl_files_dir, verbose = TRUE, ...)

# S3 method for class 'PNAAssay5'
RestorePaths(object, pxl_files_dir, verbose = TRUE, ...)

# S3 method for class 'Seurat'
RestorePaths(object, pxl_files_dir, assay = NULL, verbose = TRUE, ...)
```

## Arguments

- object:

  An object

- ...:

  Additional arguments. Currently not used.

- pxl_files_dir:

  The directory where the PXL files are stored

- verbose:

  Print messages

- assay:

  Assay name

## Value

An object with updated PXL file paths

## Details

Seurat objects created with `pixelatorR` store the MPX/PNA data in a
[`CellGraphAssay`](CellGraphAssay-class.md),
[`CellGraphAssay5`](CellGraphAssay5-class.md),
[`PNAAssay`](PNAAssay-class.md) or [`PNAAssay5`](PNAAssay5-class.md)
object. For some analytic tasks, such as component visualization, you
need to load the component graphs in memory with
[`LoadCellGraphs`](LoadCellGraphs.md).

[`LoadCellGraphs`](LoadCellGraphs.md) reads data from the edgelist(s)
stored in the PXL files(s) that the Seurat object was created from (see
[`ReadMPX_Seurat`](ReadMPX_Seurat.md) and
[`ReadPNA_Seurat`](ReadPNA_Seurat.md)). This means that the PXL files
must be accessible in the [`CellGraphAssay`](CellGraphAssay-class.md),
[`CellGraphAssay5`](CellGraphAssay5-class.md),
[`PNAAssay`](PNAAssay-class.md) or [`PNAAssay5`](PNAAssay5-class.md)
when you load the component graphs. You can check the PXL file paths
with [`FSMap`](FSMap.md).

If the PXL files are moved or copied to a different directory,
[`LoadCellGraphs`](LoadCellGraphs.md) will fail because the PXL file
paths are invalid. `RestorePaths` can be useful in these situations to
update the PXL file paths in the Seurat object. All you need is to
provide the path to the directory where the PXL files are located.

## Exported Seurat objects

A typical situation when this function is useful is when you export a
Seurat object created with `pixelatorR` to an RDS file, and share it
with someone else. To ensure that all data is available for the
recipient, you also need to share the raw PXL files and make sure that
the recipient updates the PXL file paths in the Seurat object with
`RestorePaths`.

For instance, let's assume that we have a Seurat object `seurat_obj`
created with `pixelatorR` and that the PXL files are stored in the
directory `DATA_DIR` on your local system. You can export the Seurat
object to an RDS file with:

    saveRDS(seurat_obj, file = file.path(DATA_DIR, "seurat_obj.rds"))

The content of DATA_DIR might look like this:

    DATA_DIR
    ├── seurat_obj.rds
    └── Sample01_pbmc.analysis.pxl

Now you can share the entire `DATA_DIR` folder with someone else. When
the recipient loads the RDS file in R, the PXL file paths in the Seurat
object need to be updated, because recipients `DATA_DIR` will be
different from yours. Assuming that the recipient has stored the PXL
file and the RDS object in the directory `RECIPIENT_DATA_DIR`, the
recipient can update the PXL file paths with:

    seurat_obj <- readRDS(file = file.path(RECIPIENT_DATA_DIR, "seurat_obj.rds"))
    seurat_obj <- RestorePaths(seurat_obj, pxl_files_dir = RECIPIENT_DATA_DIR)

## Examples

``` r
library(pixelatorR)
library(dplyr)

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()

# Copy PXL file to tempdir
tmp_pxl_file <- file.path(fs::path_temp(), "five_cells.pxl")
fs::file_copy(pxl_file, tmp_pxl_file)
seur_obj <- ReadMPX_Seurat(tmp_pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf45f6e47d9
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf41fa9
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf473ba3254.h5ad

# Now we can load graphs
seur_obj <- LoadCellGraphs(seur_obj, cells = colnames(seur_obj)[1])
#> ! Failed to delete temporary edge list parquet file C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf449297ec1.parquet.
#> ! Failed to delete temporary edge list parquet file C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf449297ec1.parquet.

# Removing or moving PXL file will make graphs inaccessible
fs::file_delete(tmp_pxl_file)
if (FALSE) { # \dontrun{
# This will now fail since the PXL file is missing
seur_obj <- LoadCellGraphs(seur_obj, force = TRUE)
} # }

# We can restore paths to the PXL file using RestorePaths.
# All we need is to provide the path to a directory where the PXL file is located.
pxl_files_dir <- system.file("extdata/five_cells",
  package = "pixelatorR"
)
seur_obj <- RestorePaths(seur_obj, pxl_files_dir = pxl_files_dir)

# Now LoadCellGraphs should work as expected
seur_obj <- LoadCellGraphs(seur_obj, force = TRUE)
```
