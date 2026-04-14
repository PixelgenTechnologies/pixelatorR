# Edge Rank Plot

**\[deprecated\]**

The function has been replaced by
[`MoleculeRankPlot`](MoleculeRankPlot.md) and will be removed in a
future release.

## Usage

``` r
EdgeRankPlot(object, group_by = NULL, ...)
```

## Arguments

- object:

  A Seurat object

- group_by:

  A character specifying a column to group by

- ...:

  Additional arguments to pass to
  [`MoleculeRankPlot`](MoleculeRankPlot.md)

## Examples

``` r
library(pixelatorR)

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()
seur_obj <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e86bf52238.h5ad
seur_obj
#> An object of class Seurat 
#> 80 features across 5 samples within 1 assay 
#> Active assay: mpxCells (80 features, 80 variable features)
#>  1 layer present: counts

EdgeRankPlot(seur_obj)

```
