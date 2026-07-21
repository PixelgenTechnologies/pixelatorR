# Read edgelists from a PXL file containing MPX data

This function uses arrow to read edgelists from one or several PXL
files. The edgelists are stored in parquet files in the `outdir`
directory which can be modified on disk.

## Usage

``` r
ReadMPX_arrow_edgelist(pxl_file, edge_list_file = NULL, verbose = TRUE, ...)
```

## Arguments

- pxl_file:

  Path to a PXL file

- edge_list_file:

  Path to the output edgelist.parquet file

- verbose:

  Print messages

- ...:

  Parameters passed to other methods

## Value

Nothing. The edgelist is saved to a parquet file set with
`edge_list_file`

## Examples

``` r
library(pixelatorR)

# Load example data
pxl_file <- minimal_mpx_pxl_file()
edgelist_arrow <- ReadMPX_arrow_edgelist(pxl_file)
edgelist_arrow
#> FileSystemDataset with 1 Parquet file
#> 5 columns
#> upia: string
#> upib: string
#> marker: dictionary<values=string, indices=int8>
#> count: int8
#> component: dictionary<values=string, indices=int8>
```
