# Convert polarization score table to an Assay or Assay5

Convert polarization score table to an Assay or Assay5

## Usage

``` r
PolarizationScoresToAssay(object, ...)

# S3 method for class 'data.frame'
PolarizationScoresToAssay(object, values_from = c("morans_z", "morans_i"), ...)

# S3 method for class 'MPXAssay'
PolarizationScoresToAssay(object, values_from = c("morans_z", "morans_i"), ...)

# S3 method for class 'CellGraphAssay'
PolarizationScoresToAssay(object, values_from = c("morans_z", "morans_i"), ...)

# S3 method for class 'CellGraphAssay5'
PolarizationScoresToAssay(object, values_from = c("morans_z", "morans_i"), ...)

# S3 method for class 'Seurat'
PolarizationScoresToAssay(
  object,
  assay = NULL,
  new_assay = NULL,
  values_from = c("morans_z", "morans_i"),
  ...
)
```

## Arguments

- object:

  An object with polarization scores

- ...:

  Not yet implemented

- values_from:

  What column to pick polarization scores from. Either "morans_i" or
  "morans_z"

- assay:

  Name of the [`CellGraphAssay`](CellGraphAssay-class.md) to pull
  polarization scores from

- new_assay:

  Name of the `Assay` to store the polarization scores in

## Behavior

Takes an object with polarization scores in long format and returns an
object with polarization scores in a wide format. The polarization score
table includes Moran's I and Z scores along with p-values for each
marker and component.

The wide format is an array-like object with dimensions markers x
components, where each cell is filled with a polarization score. Scores
that are missing from the polarization score table are replaced with
0's.

Different outputs are returned depending on the input object type:

- `tibble/data.frame`: returns a matrix with markers in rows and
  components in columns

- `CelGraphAssay`: returns an Assay with markers in rows and components
  in columns

- `Seurat` object: returns the Seurat object with a new Assay with
  markers in rows and components in columns

As many functions provided in Seurat works on `Assay` objects, it is
sometimes convenient to make this conversion. For instance, if we want
to compute a UMAP on the polarization scores with `RunUMAP`, we need the
values to be formatted in an `Assay`. This also makes it possible to use
various visualization functions such as `VlnPlot` or `FeaturePlor` to
show the distribution of polarization scores.

## See also

Other Spatial metrics conversion methods:
[`ColocalizationScoresToAssay()`](ColocalizationScoresToAssay.md),
[`ProximityScoresToAssay()`](ProximityScoresToAssay.md)

## Examples

``` r
library(pixelatorR)
library(SeuratObject)

# Load example data as a Seurat object
pxl_file <- minimal_mpx_pxl_file()
pol_scores <- ReadMPX_polarization(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de07423200a

# PolarizationScoresToAssay returns a matrix for a tbl_df
pol_scores_mat <- PolarizationScoresToAssay(pol_scores)
pol_scores_mat[1:4, 1:5]
#>       RCVCMP0000217 RCVCMP0000118 RCVCMP0000655 RCVCMP0000487 RCVCMP0000263
#> ACTB     -0.1335030    0.00000000    -0.2925391   -0.17816347    -0.1473378
#> B2M      -0.6759472    4.55113840     1.6401548   -3.14104176    -8.7973738
#> CD102    -0.8262732    0.03074864     0.3773469    0.03870055    -0.4982260
#> CD11a     0.1816744    0.07436839     0.8086415    0.30911581    -1.8365308

# Create a Seurat object
seur <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de03aa51e7b
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de0608d5b0b
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpInAYcT/file7de052841ef7.h5ad

# Fetch CellGraphAssay and create new polarization
# scores Assay
cg_assay <- seur[["mpxCells"]]
class(cg_assay)
#> [1] "CellGraphAssay5"
#> attr(,"package")
#> [1] "pixelatorR"
pol_assay <- PolarizationScoresToAssay(cg_assay)
class(pol_assay)
#> [1] "Assay5"
#> attr(,"package")
#> [1] "SeuratObject"

# Convert polarization scores within a Seurat object
seur <- PolarizationScoresToAssay(seur)

# After conversion, we now have a new assay called "polarization"
seur[["polarization"]]
#> Assay (v5) data with 80 features for 5 cells
#> First 10 features:
#>  ACTB, B2M, CD102, CD11a, CD11b, CD11c, CD127, CD137, CD14, CD150 
#> Layers:
#>  data 

# Switch default assay to polarization
DefaultAssay(seur) <- "polarization"

# Visualize polarization scores with Seurat
# VlnPlot(seur, features = "CD19") +
#   labs(y = "Polarization score")
```
