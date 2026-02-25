# Plot UMIs per UPIa for quality control

Plot UMIs per UPIa for quality control

## Usage

``` r
TauPlot(object, ...)

# S3 method for class 'data.frame'
TauPlot(object, group_by = NULL, ...)

# S3 method for class 'Seurat'
TauPlot(object, group_by = NULL, ...)
```

## Arguments

- object:

  A `data.frame`-like object or a `Seurat` object where columns `tau`
  and `tau_type` are present and one of the columns `umi_per_upia`,
  `mean_molecules_per_a_pixel` or `n_umi`.

- ...:

  Not yet implemented

- group_by:

  A column in the object representing a 'character' or 'factor' to group
  data by

## Value

A `ggplot` object

## See also

Other QC-plots: [`CellCountPlot()`](CellCountPlot.md)

## Examples

``` r
library(pixelatorR)

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()
seur_obj <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de021aa79c4
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de039e438b1
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpInAYcT/file7de07cab1af2.h5ad
seur_obj
#> An object of class Seurat 
#> 80 features across 5 samples within 1 assay 
#> Active assay: mpxCells (80 features, 80 variable features)
#>  1 layer present: counts

# Plot with data.frame
TauPlot(seur_obj[[]])


# Plot with Seurat object
TauPlot(seur_obj)


# Group by sample in merged data
seur_obj1 <- seur_obj2 <- seur_obj
seur_obj1$sample <- "1"
seur_obj2$sample <- "2"
seur_obj_merged <- merge(seur_obj1, seur_obj2, add.cell.ids = c("A", "B"))
TauPlot(seur_obj_merged, group_by = "sample")

```
