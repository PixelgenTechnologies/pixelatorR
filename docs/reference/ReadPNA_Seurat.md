# Load data from PNA PXL file into a `Seurat` object

This function can be used to load data from a PXL file, and returns a
`Seurat` object.

## Usage

``` r
ReadPNA_Seurat(
  pxl_file,
  assay = "PNA",
  return_pna_assay = TRUE,
  load_proximity_scores = TRUE,
  calc_log2_ratio = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- pxl_file:

  Path to a PXL file

- assay:

  Assay name

- return_pna_assay:

  Logical specifying whether the count data should be stored in a
  `PNAAssay`/`PNAAssay5`. If set to `FALSE`, the count data will be
  stored in a `Assay`/`Assay5` instead. For the latter case, many
  features provided in `pixelatorR` will be unavailable. This can be
  useful if you only intend to analyze the abundance data.

- load_proximity_scores:

  Logical specifying whether the proximity scores should be loaded into
  the `PNAAssay`/`PNAAssay5`. If you only intend to analyze abundance
  data or PNA graphs, you can set this parameter to `FALSE` to use less
  memory. This parameter only have an effect if
  `return_pna_assay = TRUE`.

- calc_log2_ratio:

  A logical specifying whether to calculate and add a log2ratio column
  to the output table. Default is `TRUE`

- verbose:

  Print messages

- ...:

  Additional parameters passed to
  [`CreateSeuratObject`](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html)

## Value

An object of class `Seurat`

## Details

By default, the count matrix is returned in a `PNAAssay` object. If the
global option `Seurat.object.assay.version` is set to `"v5"`, the
function will return a `PNAAssay5` object instead.

## See also

Other PXL-data-loaders: [`ReadPNA_counts()`](ReadPNA_counts.md),
[`ReadPNA_metadata()`](ReadPNA_metadata.md)

## Examples

``` r
library(pixelatorR)

# Crete example Seurat object
pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
seur_obj
#> An object of class Seurat 
#> 158 features across 5 samples within 1 assay 
#> Active assay: PNA (158 features, 158 variable features)
#>  1 layer present: counts
```
