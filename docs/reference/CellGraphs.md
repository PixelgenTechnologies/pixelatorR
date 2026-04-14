# CellGraphs

Get and set [`CellGraph`](CellGraph-class.md) lists for different
objects.

## Usage

``` r
CellGraphs(object, ...)

CellGraphs(object, ...) <- value

# S3 method for class 'MPXAssay'
CellGraphs(object, ...)

# S3 method for class 'MPXAssay'
CellGraphs(object, ...) <- value

# S3 method for class 'PNAAssay'
CellGraphs(object, ...)

# S3 method for class 'PNAAssay5'
CellGraphs(object, ...)

# S3 method for class 'PNAAssay'
CellGraphs(object, ...) <- value

# S3 method for class 'PNAAssay5'
CellGraphs(object, ...) <- value

# S3 method for class 'Seurat'
CellGraphs(object, ...)

# S3 method for class 'Seurat'
CellGraphs(object, ...) <- value
```

## Arguments

- object:

  An object with cellgraphs

- ...:

  Additional arguments

- value:

  A named list with [`CellGraph`](CellGraph-class.md) objects to replace
  the current cellgraphs

## Value

Returns a list of [`CellGraph`](CellGraph-class.md) objects. If there
are no [`CellGraph`](CellGraph-class.md) objects present, returns a
named list where each element is `NULL`.

## See also

[`PolarizationScores()`](PolarizationScores.md) and
[`ColocalizationScores()`](ColocalizationScores.md) for getting/setting
spatial metrics

## Examples

``` r
library(pixelatorR)
library(dplyr)
library(tidygraph)

pxl_file <- minimal_mpx_pxl_file()
counts <- ReadMPX_counts(pxl_file)
#> ℹ Loading count data from C:/Users/max/AppData/Local/R/win-library/4.5/pixelatorR/extdata/five_cells/five_cells.pxl
edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#> ℹ Loading item(s) from: C:/Users/max/AppData/Local/R/win-library/4.5/pixelatorR/extdata/five_cells/five_cells.pxl
#> →   Loading edgelist data
#> ✔ Returning a 'tbl_df' object
components <- colnames(counts)
edgelist_split <-
  edgelist %>%
  select(upia, upib, component) %>%
  distinct() %>%
  group_by(component) %>%
  group_split() %>%
  setNames(nm = components)

# Convert data into a list of CellGraph objects
bipartite_graphs <- lapply(edgelist_split, function(x) {
  x <- x %>% as_tbl_graph(directed = FALSE)
  x <- x %>% mutate(node_type = case_when(name %in% edgelist$upia ~ "A", TRUE ~ "B"))
  attr(x, "type") <- "bipartite"
  CreateCellGraphObject(cellgraph = x)
})

# CellGraphs getter CellGraphAssay
# ---------------------------------

# Create CellGraphAssay
cg_assay <- CreateCellGraphAssay(counts = counts, cellgraphs = bipartite_graphs)
cg_assay
#> CellGraphAssay data with 80 features for 5 cells
#> First 10 features:
#>  CD274, CD44, CD25, CD279, CD41, HLA-ABC, CD54, CD26, CD27, CD38 
#> Loaded CellGraph objects:
#>  5

# Get cellgraphs from a CellGraphAssay object
CellGraphs(cg_assay)
#> $RCVCMP0000217
#> A CellGraph object containing a bipartite graph with 3507 nodes and 7580 edges
#> 
#> $RCVCMP0000118
#> A CellGraph object containing a bipartite graph with 2470 nodes and 5138 edges
#> 
#> $RCVCMP0000487
#> A CellGraph object containing a bipartite graph with 4225 nodes and 8150 edges
#> 
#> $RCVCMP0000655
#> A CellGraph object containing a bipartite graph with 4340 nodes and 9918 edges
#> 
#> $RCVCMP0000263
#> A CellGraph object containing a bipartite graph with 2258 nodes and 4303 edges
#> 


# CellGraphs setter CellGraphAssay
# ---------------------------------

# Set cellgraphs in a CellGraphAssay object
CellGraphs(cg_assay) <- cg_assay@cellgraphs

library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
CellGraphs(seur_obj[["PNA"]])
#> $`0a45497c6bfbfb22`
#> NULL
#> 
#> $`2708240b908e2eba`
#> NULL
#> 
#> $c3c393e9a17c1981
#> NULL
#> 
#> $d4074c845bb62800
#> NULL
#> 
#> $efe0ed189cb499fc
#> NULL
#> 

# Set cellgraphs in a PNAAssay object
CellGraphs(seur_obj[["PNA"]]) <- CellGraphs(seur_obj[["PNA"]])


# CellGraphs getter Seurat
# ---------------------------------
pxl_file <- minimal_mpx_pxl_file()
se <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpmOhqBt/file62e86d6d4f07.h5ad

# Get cellgraphs from a Seurat object
CellGraphs(se)
#> $RCVCMP0000217
#> NULL
#> 
#> $RCVCMP0000118
#> NULL
#> 
#> $RCVCMP0000487
#> NULL
#> 
#> $RCVCMP0000655
#> NULL
#> 
#> $RCVCMP0000263
#> NULL
#> 

# CellGraphs setter Seurat
# ---------------------------------

# Set cellgraphs in a Seurat object
CellGraphs(se) <- cg_assay@cellgraphs
```
