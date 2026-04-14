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

  Information to pull from object (cellgraph, meta_data, layout)

- value:

  A new variable to place in `slot`

## Value

`GetCellGraphData`: A [`CellGraph`](CellGraph-class.md) object slot

`CellGraphData<-`: A [`CellGraph`](CellGraph-class.md) with updated data

## Examples

``` r
library(pixelatorR)
library(dplyr)
library(tidygraph)

edge_list <-
  ReadMPX_item(
    minimal_mpx_pxl_file(),
    items = "edgelist"
  )
#> ℹ Loading item(s) from: C:/Users/max/AppData/Local/R/win-library/4.5/pixelatorR/extdata/five_cells/five_cells.pxl
#> →   Loading edgelist data
#> ✔ Returning a 'tbl_df' object
bipart_graph <-
  edge_list %>%
  select(upia, upib, marker) %>%
  distinct() %>%
  as_tbl_graph(directed = FALSE) %>%
  mutate(node_type = case_when(name %in% edge_list$upia ~ "A", TRUE ~ "B"))
attr(bipart_graph, "type") <- "bipartite"

cg <- CreateCellGraphObject(cellgraph = bipart_graph)

# Get slot data
CellGraphData(cg, slot = "cellgraph")
#> # A tbl_graph: 16800 nodes and 68255 edges
#> #
#> # An undirected multigraph with 5 components
#> #
#> # Node Data: 16,800 × 2 (active)
#>    name                      node_type
#>    <chr>                     <chr>    
#>  1 GGTATATATTTTAAGTAGTTTAGTA A        
#>  2 ATTATGTTGTTGTATGTTTATTATT A        
#>  3 TCCGTGGTTTGATACTCGGGAATTT A        
#>  4 AGTGTAAGAGGTTGTTTCTTAGAAA A        
#>  5 GAGCAGACAATGGCGCTTAGCTAAA A        
#>  6 CAATTTTGACCTAGTTGTGGCCAAG A        
#>  7 ATGTCAGAGTGAATAGTTGTGATTG A        
#>  8 ATAAAGTGACCAGTTGAATGGGCCC A        
#>  9 TGCCTCTTGAGCTGTTGATCGCCTG A        
#> 10 ATTAGATGTAGAGCGCGAAAATTGG A        
#> # ℹ 16,790 more rows
#> #
#> # Edge Data: 68,255 × 3
#>    from    to marker
#>   <int> <int> <fct> 
#> 1     1 11034 CD27  
#> 2     2 11035 CD27  
#> 3     3 11036 CD27  
#> # ℹ 68,252 more rows

# Set slot data
CellGraphData(cg, slot = "cellgraph") <- CellGraphData(cg, slot = "cellgraph")
```
