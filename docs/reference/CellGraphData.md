# Get and set CellGraph object data

Get and set CellGraph object data

## Usage

``` r
CellGraphData(object, slot = "cellgraph")

CellGraphData(object, slot = "cellgraph") <- value
```

## Arguments

- object:

  A [`CellGraph`](CellGraph-class.md) object

- slot:

  Information to pull from object (`cellgraph`, `counts`, `layout`,
  `layers`, `meta.data`, `reductions`). `meta_data` is accepted as an
  alias for `meta.data`.

- value:

  A new variable to place in `slot`

## Value

`GetCellGraphData`: A [`CellGraph`](CellGraph-class.md) object slot

`CellGraphData<-`: A [`CellGraph`](CellGraph-class.md) with updated data

## Examples

``` r
library(pixelatorR)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
library(tidygraph)
#> 
#> Attaching package: ‘tidygraph’
#> The following object is masked from ‘package:stats’:
#> 
#>     filter

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
cg <- CellGraphs(se)[[1]]

# Get slot data
CellGraphData(cg, slot = "cellgraph")
#> # A tbl_graph: 43543 nodes and 97014 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 43,543 × 2 (active)
#>    name                   node_type
#>    <chr>                  <chr>    
#>  1 61208583141770358-umi1 umi1     
#>  2 50526950249468550-umi2 umi2     
#>  3 69733109123764664-umi1 umi1     
#>  4 43235234960499656-umi2 umi2     
#>  5 16002757515879905-umi1 umi1     
#>  6 4606209975865882-umi2  umi2     
#>  7 59822389138925142-umi1 umi1     
#>  8 28470384711643019-umi2 umi2     
#>  9 64270251753030037-umi1 umi1     
#> 10 63550490663717685-umi2 umi2     
#> # ℹ 43,533 more rows
#> #
#> # Edge Data: 97,014 × 2
#>    from    to
#>   <int> <int>
#> 1     1     2
#> 2     3     4
#> 3     5     6
#> # ℹ 97,011 more rows

# Set slot data
CellGraphData(cg, slot = "cellgraph") <- CellGraphData(cg, slot = "cellgraph")
```
