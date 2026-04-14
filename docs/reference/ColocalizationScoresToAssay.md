# Convert colocalization score table to an Assay or Assay5

Convert colocalization score table to an Assay or Assay5

## Usage

``` r
ColocalizationScoresToAssay(object, ...)

# S3 method for class 'data.frame'
ColocalizationScoresToAssay(
  object,
  values_from = c("pearson_z", "pearson"),
  ...
)

# S3 method for class 'MPXAssay'
ColocalizationScoresToAssay(
  object,
  values_from = c("pearson_z", "pearson"),
  ...
)

# S3 method for class 'CellGraphAssay'
ColocalizationScoresToAssay(
  object,
  values_from = c("pearson_z", "pearson"),
  ...
)

# S3 method for class 'CellGraphAssay5'
ColocalizationScoresToAssay(
  object,
  values_from = c("pearson_z", "pearson"),
  ...
)

# S3 method for class 'Seurat'
ColocalizationScoresToAssay(
  object,
  assay = NULL,
  new_assay = NULL,
  values_from = c("pearson_z", "pearson"),
  ...
)
```

## Arguments

- object:

  An object with colocalization scores

- ...:

  Not yet implemented

- values_from:

  What column to pick colocalization scores from. One of "pearson" or
  "pearson_z"

- assay:

  Name of the [`CellGraphAssay`](CellGraphAssay-class.md) to pull
  polarization scores from

- new_assay:

  Name of the `Assay` to store the polarization scores in

## Behavior

Takes an object with colocalization scores in long format and returns an
object with colocalization scores in a wide format. The colocalization
score table includes various colocalization scores along with p-values
for each pair markers and component.

The wide format is an array-like object with dimensions (markers_1 \*
marker_2) x components, where each cell is filled with a polarization
score. Scores that are missing from the colocalization score table are
replaced with 0's.

Different outputs are returned depending on the input object type:

- `tibble/data.frame`: returns a matrix with marker pairs in rows and
  components in columns

- `CelGraphAssay`: returns an Assay with marker pairs in rows and
  components in columns

- `Seurat` object: returns the Seurat object with a new Assay with
  marker pairs in rows and components in columns

As many functions provided in Seurat works on `Assay` objects, it is
sometimes convenient to make this conversion. For instance, if we want
to compute a UMAP on the colocalization scores with `RunUMAP`, we need
the values to be formatted in an `Assay`. This also makes it possible to
use various visualization functions such as `VlnPlot` or `FeaturePlor`
to show the distribution of colocalization scores.

## See also

Other Spatial metrics conversion methods:
[`PolarizationScoresToAssay()`](PolarizationScoresToAssay.md),
[`ProximityScoresToAssay()`](ProximityScoresToAssay.md)

## Examples

``` r
library(pixelatorR)
library(SeuratObject)

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()
col_scores <- ReadMPX_colocalization(pxl_file)
#> ℹ Loading item(s) from: C:/Users/max/AppData/Local/R/win-library/4.5/pixelatorR/extdata/five_cells/five_cells.pxl
#> →   Loading colocalization data
#> ✔ Returning a 'tbl_df' object

# ColocalizationScoresToAssay returns a matrix for a tbl_df
col_scores_mat <- ColocalizationScoresToAssay(col_scores)
col_scores_mat[1:4, 1:5]
#>            RCVCMP0000217 RCVCMP0000118 RCVCMP0000487 RCVCMP0000655
#> ACTB/B2M      -0.9401321      0.000000     -2.879992     4.5740882
#> ACTB/CD102     1.9241318      0.000000      3.157566    -0.4834448
#> B2M/CD102      4.7268224      2.190192     -3.101690     2.3802011
#> ACTB/CD11a     2.1011979      0.000000     -1.804651     8.0947966
#>            RCVCMP0000263
#> ACTB/B2M       1.9368818
#> ACTB/CD102     0.5934362
#> B2M/CD102     -5.2602955
#> ACTB/CD11a     5.4727656

# Create a Seurat object
seur <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e81355d37.h5ad

# Fetch CellGraphAssay and create new polarization
# scores Assay
cg_assay <- seur[["mpxCells"]]
class(cg_assay)
#> [1] "CellGraphAssay5"
#> attr(,"package")
#> [1] "pixelatorR"
col_assay <- ColocalizationScoresToAssay(cg_assay)
class(col_assay)
#> [1] "Assay5"
#> attr(,"package")
#> [1] "SeuratObject"

# Convert colocalization scores within a Seurat object
seur <- ColocalizationScoresToAssay(seur)

# After conversion, we now have a new assay called "colocalization"
seur[["colocalization"]]
#> Assay (v5) data with 3160 features for 5 cells
#> First 10 features:
#>  ACTB/B2M, ACTB/CD102, B2M/CD102, ACTB/CD11a, B2M/CD11a, CD102/CD11a,
#> ACTB/CD11b, B2M/CD11b, CD102/CD11b, CD11a/CD11b 
#> Layers:
#>  data 

# Switch default assay to polarization
DefaultAssay(seur) <- "colocalization"

# Visualize colocalization scores with Seurat
# VlnPlot(seur, features = "CD19") +
#   ggplot2::labs(y = "Colocalization score")
```
