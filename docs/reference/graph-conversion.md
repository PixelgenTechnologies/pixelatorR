# A-node projection

**\[deprecated\]**

A function to create an A node projection graph from an edgelist
representing a bipartite graph.

## Usage

``` r
edgelist_to_simple_Anode_graph(object, ...)

# S3 method for class 'data.frame'
edgelist_to_simple_Anode_graph(
  object,
  components = NULL,
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'FileSystemDataset'
edgelist_to_simple_Anode_graph(object, components = NULL, verbose = TRUE, ...)
```

## Arguments

- object:

  An object containing an edgelist

- ...:

  Not yet implemented

- components:

  Components to compute A node projection for

- cl:

  A cluster object created by makeCluster, or an integer to indicate
  number of child-processes (integer values are ignored on Windows) for
  parallel evaluations. See Details on performance in the documentation
  for `pbapply`. The default is NULL, which means that no
  parallelization is used.

- verbose:

  Print messages

## Value

An A-node-projected graph

## Examples

``` r
library(pixelatorR)
library(tibble)
#> Warning: package 'tibble' was built under R version 4.5.3
#> 
#> Attaching package: 'tibble'
#> The following object is masked from 'package:igraph':
#> 
#>     as_data_frame

pxl_file <- minimal_mpx_pxl_file()

# Load edgelist
el <- ReadMPX_arrow_edgelist(pxl_file)
#> ℹ Extracting edgelist.parquet file to C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/edgelist.parquet
#> ✔ Returning FileSystemDataset

# Convert to tbl_df
el_tbl_df <- as_tibble(el)

# data.frame method --------------------------------
# Create list of A-node projected tbl_graphs
anode_proj <- edgelist_to_simple_Anode_graph(el_tbl_df)
#> ℹ Simplifying edge list
#> ℹ Creating A-node projected graphs
#> ✔ Returning an A-node projected graphs

# FileSystemDataset method -------------------------
# Create list of A-node projected tbl_graphs
anode_proj <- edgelist_to_simple_Anode_graph(el)
#> ℹ Simplifying edge list
#> ℹ Creating A-node projected graphs
#> ✔ Returning an A-node projected graphs
```
