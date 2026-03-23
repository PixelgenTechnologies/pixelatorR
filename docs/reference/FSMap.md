# FS map

Get and set `fs_map` tibble for a
[`CellGraphAssay`](CellGraphAssay-class.md),
[`CellGraphAssay5`](CellGraphAssay5-class.md) or a `Seurat` object

## Usage

``` r
FSMap(object, ...)

FSMap(object, ...) <- value

# S3 method for class 'MPXAssay'
FSMap(object, ...)

# S3 method for class 'MPXAssay'
FSMap(object, ...) <- value

# S3 method for class 'PNAAssay'
FSMap(object, ...)

# S3 method for class 'PNAAssay5'
FSMap(object, ...)

# S3 method for class 'PNAAssay'
FSMap(object, ...) <- value

# S3 method for class 'PNAAssay5'
FSMap(object, ...) <- value

# S3 method for class 'Seurat'
FSMap(object, ...)

# S3 method for class 'Seurat'
FSMap(object, ...) <- value
```

## Arguments

- object:

  A `PNAAssay` or `PNAAssay5` object

- ...:

  Additional arguments (not used)

- value:

  A new `tbl_df` to replace the current fs_map with

## Value

`FSMap`: An `fs_map` tibble

`FSMap<-`: An object with an updated `fs_map` tibble

## See also

[`CellGraphs()`](CellGraphs.md) for getting/setting
[`CellGraph`](CellGraph-class.md) lists and
[`PolarizationScores()`](PolarizationScores.md),[`ColocalizationScores()`](ColocalizationScores.md)
for getting/setting spatial metrics

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
FSMap(seur_obj[["PNA"]])
#> # A tibble: 1 × 3
#>   id_map           sample pxl_file                                              
#>   <list>            <int> <chr>                                                 
#> 1 <tibble [5 × 2]>      1 "C:\\Users\\max\\AppData\\Local\\Temp\\Rtmp2f7Zq5\\te…

# If the PXL has been moved, we can update the fs_map
# Here we copy the test PXL file to a temporary location
# to illustrate how to update the fs_map
temp_file <- fs::file_temp(ext = "pxl")
fs::file_copy(path = pxl_file, temp_file)

# Change file path in fs_map
FSMap(seur_obj[["PNA"]]) <- FSMap(seur_obj[["PNA"]]) %>%
  mutate(pxl_file = temp_file)

# Now the fs_map has been updated with the correct path
FSMap(seur_obj[["PNA"]])
#> # A tibble: 1 × 3
#>   id_map           sample pxl_file                                              
#>   <list>            <int> <fs::path>                                            
#> 1 <tibble [5 × 2]>      1 …ax/AppData/Local/Temp/Rtmpampkmn/file5bf472535cc9.pxl

pxl_file <- minimal_mpx_pxl_file()
seur_obj <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf43cf641bb
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf419136800
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf43131138c.h5ad

# Check PXL file paths in a Seurat object
FSMap(seur_obj)
#> # A tibble: 1 × 3
#>   id_map           sample pxl_file                                              
#>   <list>            <int> <chr>                                                 
#> 1 <tibble [5 × 2]>      1 "C:\\Users\\max\\AppData\\Local\\Temp\\Rtmp2f7Zq5\\te…

library(pixelatorR)

# Create example Seurat object
pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins

# Replace FSMap in Seurat object
FSMap(seur_obj) <- FSMap(seur_obj) %>%
  mutate(pxl_file = fs::path_rel(pxl_file))
```
