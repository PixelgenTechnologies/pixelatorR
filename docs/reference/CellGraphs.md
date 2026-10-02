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

Returns a [`CellGraphList`](CellGraphList.md). Unloaded graphs are
stored as `NULL` elements.

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
#> ℹ Loading count data from /tmp/Rtmpoj0cDG/temp_libpath4d102621131d/pixelatorR/extdata/five_cells/five_cells.pxl
edgelist <- ReadMPX_item(pxl_file, items = "edgelist")
#> ℹ Loading item(s) from: /tmp/Rtmpoj0cDG/temp_libpath4d102621131d/pixelatorR/extdata/five_cells/five_cells.pxl
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
#> A CellGraphList with 5 loaded CellGraph object(s) out of 5 
#> Names: RCVCMP0000217, RCVCMP0000118, RCVCMP0000487, RCVCMP0000655, RCVCMP0000263  


# CellGraphs setter CellGraphAssay
# ---------------------------------

# Set cellgraphs in a CellGraphAssay object
CellGraphs(cg_assay) <- cg_assay@cellgraphs

library(pixelatorR)

pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpjKKHFf/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
CellGraphs(seur_obj[["PNA"]])
#> A CellGraphList with 0 loaded CellGraph object(s) out of 5 
#> Names: 0a45497c6bfbfb22, 2708240b908e2eba, c3c393e9a17c1981, d4074c845bb62800, efe0ed189cb499fc  

# Set cellgraphs in a PNAAssay object
CellGraphs(seur_obj[["PNA"]]) <- CellGraphs(seur_obj[["PNA"]])


# CellGraphs getter Seurat
# ---------------------------------
pxl_file <- minimal_mpx_pxl_file()
se <- ReadMPX_Seurat(pxl_file)
#> ✔ Created a 'Seurat' object with 5 cells and 80 targeted surface proteins

# Get cellgraphs from a Seurat object
CellGraphs(se)
#> A CellGraphList with 0 loaded CellGraph object(s) out of 5 
#> Names: RCVCMP0000217, RCVCMP0000118, RCVCMP0000487, RCVCMP0000655, RCVCMP0000263  

# CellGraphs setter Seurat
# ---------------------------------

# Set cellgraphs in a Seurat object
CellGraphs(se) <- cg_assay@cellgraphs
```
