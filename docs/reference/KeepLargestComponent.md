# Keep largest component

Finds connected components of a graph and returns the largest component

## Usage

``` r
KeepLargestComponent(object, ...)

# S3 method for class 'tbl_graph'
KeepLargestComponent(object, verbose = TRUE, ...)

# S3 method for class 'CellGraph'
KeepLargestComponent(object, verbose = TRUE, ...)

# S3 method for class 'MPXAssay'
KeepLargestComponent(object, verbose = TRUE, ...)

# S3 method for class 'PNAAssay'
KeepLargestComponent(object, verbose = TRUE, ...)

# S3 method for class 'PNAAssay5'
KeepLargestComponent(object, verbose = TRUE, ...)

# S3 method for class 'Seurat'
KeepLargestComponent(object, assay = NULL, verbose = TRUE, ...)
```

## Arguments

- object:

  An object

- ...:

  Parameters passed to other methods

- verbose:

  Print messages

- assay:

  Name of a `CellGraphAssay`, `CellGraphAssay5`, `PNAAssay`
  or`PNAAssay5` assay stored on the `Seurat` object

## Examples

``` r
library(pixelatorR)
library(tidygraph)

pxl_file <- minimal_mpx_pxl_file()

# Read edgelist
edgelist <- ReadMPX_arrow_edgelist(pxl_file)
#> ℹ Extracting edgelist.parquet file to C:/Users/max/AppData/Local/Temp/RtmpInAYcT/edgelist.parquet
#> ✔ Returning FileSystemDataset

# Load graph from edge list and store in a CellGraph object
cg <- LoadCellGraphs(edgelist, cells = "RCVCMP0000217", data_type = "MPX")[[1]]
#> →   Loading 1 edgelist(s) as bipartite graph(s)
cg
#> A CellGraph object containing a bipartite graph with 2470 nodes and 5138 edges
#> Number of markers:  79 

# Fetch tbl_graph from CellGraph object
g <- CellGraphData(cg, slot = "cellgraph")
g
#> # A tbl_graph: 2470 nodes and 5138 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 2,470 × 2 (active)
#>    name                        node_type
#>    <chr>                       <chr>    
#>  1 GGTATATATTTTAAGTAGTTTAGTA-A A        
#>  2 ATTATGTTGTTGTATGTTTATTATT-A A        
#>  3 TCCGTGGTTTGATACTCGGGAATTT-A A        
#>  4 AGTGTAAGAGGTTGTTTCTTAGAAA-A A        
#>  5 GAGCAGACAATGGCGCTTAGCTAAA-A A        
#>  6 CAATTTTGACCTAGTTGTGGCCAAG-A A        
#>  7 ATGTCAGAGTGAATAGTTGTGATTG-A A        
#>  8 ATAAAGTGACCAGTTGAATGGGCCC-A A        
#>  9 TGCCTCTTGAGCTGTTGATCGCCTG-A A        
#> 10 ATTAGATGTAGAGCGCGAAAATTGG-A A        
#> # ℹ 2,460 more rows
#> #
#> # Edge Data: 5,138 × 2
#>    from    to
#>   <int> <int>
#> 1     1  1396
#> 2     1  1599
#> 3     1  1734
#> # ℹ 5,135 more rows

# Break graph by removing random edges
set.seed(132)
g <- g %E>%
  filter(from %in% sample(from, n() - 500))

# Fetch largest component from a tbl_graph
g_largest <- KeepLargestComponent(g)
#> ℹ Removed 35 out of 2470 nodes
g_largest
#> # A tbl_graph: 2435 nodes and 5103 edges
#> #
#> # An undirected simple graph with 1 component
#> #
#> # Node Data: 2,435 × 2 (active)
#>    name                        node_type
#>    <chr>                       <chr>    
#>  1 GGTATATATTTTAAGTAGTTTAGTA-A A        
#>  2 ATTATGTTGTTGTATGTTTATTATT-A A        
#>  3 TCCGTGGTTTGATACTCGGGAATTT-A A        
#>  4 AGTGTAAGAGGTTGTTTCTTAGAAA-A A        
#>  5 GAGCAGACAATGGCGCTTAGCTAAA-A A        
#>  6 CAATTTTGACCTAGTTGTGGCCAAG-A A        
#>  7 ATGTCAGAGTGAATAGTTGTGATTG-A A        
#>  8 ATAAAGTGACCAGTTGAATGGGCCC-A A        
#>  9 TGCCTCTTGAGCTGTTGATCGCCTG-A A        
#> 10 ATTAGATGTAGAGCGCGAAAATTGG-A A        
#> # ℹ 2,425 more rows
#> #
#> # Edge Data: 5,103 × 2
#>    from    to
#>   <int> <int>
#> 1     1  1362
#> 2     1  1565
#> 3     1  1700
#> # ℹ 5,100 more rows

# Add new graph to CellGraph object
cg@cellgraph <- g_largest

# Fetch largest component from a CellGraph
cg_largest <- KeepLargestComponent(cg)
#> ℹ Graph is already connected
cg_largest
#> A CellGraph object containing a bipartite graph with 2435 nodes and 5103 edges
#> Number of markers:  79 
```
