# Load data from PXL file into a `Seurat` object

This wrapper function can be used to load data from a PXL file, and
returns a `Seurat` object.

## Usage

``` r
ReadMPX_Seurat(
  filename,
  assay = "mpxCells",
  return_cellgraphassay = TRUE,
  load_cell_graphs = FALSE,
  load_polarity_scores = TRUE,
  load_colocalization_scores = TRUE,
  add_additional_assays = FALSE,
  edgelist_outdir = NULL,
  overwrite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- filename:

  Path to a PXL file

- assay:

  Assay name

- return_cellgraphassay:

  Should data be loaded as a `CellGraphAssay` object?

- load_cell_graphs:

  Should the cellgraphs be loaded into the `CellGraphAssay` object?

- load_polarity_scores, load_colocalization_scores:

  Logical specifying if the polarity and colocalization scores should be
  loaded. These parameters only have an effect if
  `return_cellgraphassay = TRUE`.

- add_additional_assays:

  If other matrix representations are stored in the PXL file, for
  instance CLR-normalized counts or denoised, set this parameter to
  `TRUE` to load these in separate Assays.

- edgelist_outdir:

  A directory where the edgelist should be stored

- overwrite:

  Should `edgelist_outdir` be overwritten?

- verbose:

  Print messages

- ...:

  Additional parameters passed to
  [`CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html)

## Value

An object of class `Seurat`

## Details

By default, the MPX count matrix is returned in a `CellGraphAssay`
object. Graphs are not loaded directly unless `load_cell_graphs = TRUE`.
Graphs can also be loaded at a later stage with
[`LoadCellGraphs`](LoadCellGraphs.md).

When setting the global option `Seurat.object.assay.version` to `"v5"`,
the function will return a `CellGraphAssay5` object instead.

## See also

Other data-loaders: [`ReadMPX_item()`](ReadMPX_item.md),
[`ReadPNA_proximity()`](ReadPNA_proximity.md)

## Examples

``` r
library(pixelatorR)

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()
seur_obj <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de04dbb17e1
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de035fd8f5
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpInAYcT/file7de065317336.h5ad
seur_obj
#> An object of class Seurat 
#> 80 features across 5 samples within 1 assay 
#> Active assay: mpxCells (80 features, 80 variable features)
#>  1 layer present: counts
```
