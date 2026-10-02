# CellGraph Methods

Methods for [`CellGraph`](CellGraph-class.md) objects for generics
defined in other packages

## Usage

``` r
# S3 method for class 'CellGraph'
Layers(object, search = NULL, ...)

# S3 method for class 'CellGraph'
LayerData(object, layer = "counts", ...)

# S3 method for class 'CellGraph'
LayerData(object, layer = "counts", ...) <- value

# S3 method for class 'CellGraph'
Embeddings(object, reduction = NULL, ...)

# S3 method for class 'CellGraph'
Loadings(object, reduction = NULL, ...)

# S3 method for class 'CellGraph'
Stdev(object, reduction = NULL, ...)

# S3 method for class 'CellGraph'
AddMetaData(object, metadata, col.name = NULL, ...)

# S3 method for class 'CellGraph'
FetchData(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  clean = TRUE,
  add_marker = FALSE,
  ...
)

# S4 method for class 'CellGraph'
show(object)

# S4 method for class 'CellGraph,character,missing'
x[[i, j, ..., drop = TRUE]]

# S4 method for class 'CellGraph,character,missing,ANY'
x[[i, j, ...]] <- value

# S3 method for class 'CellGraph'
subset(x, nodes, ...)

# S3 method for class 'CellGraphList'
FetchData(
  object,
  vars,
  cells = NULL,
  layer = NULL,
  clean = FALSE,
  add_marker = FALSE,
  ...
)
```

## Arguments

- object:

  A [`CellGraph`](CellGraph-class.md) or
  [`CellGraphList`](CellGraphList.md) object

- search:

  Optional layer name or pattern passed to
  [`Layers`](https://satijalab.github.io/seurat-object/reference/Layers.html)

- ...:

  Currently not used

- layer:

  Name of a node matrix layer. Use `"counts"` for the counts slot. For
  `FetchData`, `NULL` (default) selects `"counts"` when present,
  otherwise the first extra layer.

- value:

  Replacement value

- reduction:

  Name of a stored [`NodeDimReduc`](NodeDimReduc-class.md). Defaults to
  the first reduction when `NULL`.

- metadata:

  A vector, matrix, or `data.frame` of node metadata. Nodes are matched
  by name, so metadata may cover a subset of the graph; the remaining
  nodes get `NA`. Names that are not graph nodes are dropped.

- col.name:

  Name of the metadata column when `metadata` is a vector

- vars:

  Variables to fetch: marker names, node metadata columns, graph vertex
  attributes, or reduction embedding columns (for example `"PC_1"`). A
  variable missing from the graph (or from every cell in a list) aborts.
  A variable present on some cells in a list and missing on others is
  filled with `NA` and a warning names those cells.

- cells:

  For `FetchData.CellGraph`, nodes to collect (default is all nodes).
  Numeric indices are allowed, matching
  [`FetchData`](https://satijalab.github.io/seurat-object/reference/FetchData.html).
  For `FetchData.CellGraphList`, component IDs (default is all loaded
  graphs). Unloaded graphs raise an error when they are included in
  `cells`.

- clean:

  If `TRUE`, remove nodes that are missing data for every requested
  variable. `FetchData.CellGraph` defaults to `TRUE`.
  `FetchData.CellGraphList` defaults to `FALSE` so graphs that lack the
  requested variables still appear with `NA` values. A `marker` column
  added with `add_marker` is not treated as a requested variable.

- add_marker:

  If `TRUE`, add a `marker` column with the marker label of each node
  from the one-hot counts matrix. Nodes with no count are `NA`. `vars`
  cannot include `marker` when this is `TRUE`.

- x:

  A [`CellGraph`](CellGraph-class.md) object

- i:

  Name of a stored reduction

- j, drop:

  Required by the S4 `[[` generic and ignored

- nodes:

  A character vector of node names

## Value

`FetchData.CellGraph`: a `data.frame` with nodes as rows and requested
variables as columns. `FetchData.CellGraphList`: a `data.frame` with a
`component` column identifying the source graph and the requested
variables. Row names are node IDs and must be unique across the combined
graphs. `add_marker = TRUE` also adds a `marker` column with the marker
label of each node. `subset`: a `CellGraph` object containing only the
specified nodes.

## Details

Variable names must be unique across the graph node table, `meta.data`,
reduction embeddings, and matrix features. Count and layer matrices may
share feature names because `layer` selects the matrix to search.

## Functions

- `FetchData(CellGraph)`: Pull node-level data from a `CellGraph`

- `show(CellGraph)`: Show a `CellGraph` object

- `x[[i`: Extract a `NodeDimReduc` by name

- `` `[[`(x = CellGraph, i = character, j = missing) <- value ``: Add or
  replace a `NodeDimReduc`

- `subset(CellGraph)`: Subset a `CellGraph` object

- `FetchData(CellGraphList)`: Pull node-level data from each loaded
  `CellGraph` in a `CellGraphList`. Unlike
  [`FetchLayoutData`](FetchLayoutData.md), this does not require a
  stored layout and does not reserve coordinate names, so `vars` may
  include `x`, `y`, or `z` when those columns exist on the graphs.
  `component` is reserved for the source graph ID. Node IDs are used as
  row names and must be unique across the graphs being combined.
  Variables missing from a graph are filled with `NA` and a warning
  names the cells that lacked the variable. Variables missing from every
  graph abort. Graphs that do not have a requested `layer` omit only
  features from that layer; metadata, vertex attributes, and reductions
  are kept. Remaining names are not looked up in `counts` or other
  layers. `clean` defaults to `FALSE` so those missing values are kept.
  `add_marker = TRUE` adds a `marker` column from the one-hot counts
  matrix of each graph.

## Examples

``` r
library(pixelatorR)
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

# Show method
cg
#> A CellGraph object containing a bipartite graph with 43543 nodes and 97014 edges
#> Number of markers:  149 

# Fetch marker counts, node attributes, or embeddings
head(SeuratObject::FetchData(cg, vars = colnames(cg@counts)[1]))
#>                        B2M
#> 61208583141770358-umi1   0
#> 50526950249468550-umi2   0
#> 69733109123764664-umi1   0
#> 43235234960499656-umi2   0
#> 16002757515879905-umi1   0
#> 4606209975865882-umi2    0

# FetchData.CellGraphList combines loaded graphs without requiring a layout
cgl <- CellGraphs(se)
head(SeuratObject::FetchData(cgl, vars = colnames(cg@counts)[1]))
#>                               component B2M
#> 61208583141770358-umi1 0a45497c6bfbfb22   0
#> 50526950249468550-umi2 0a45497c6bfbfb22   0
#> 69733109123764664-umi1 0a45497c6bfbfb22   0
#> 43235234960499656-umi2 0a45497c6bfbfb22   0
#> 16002757515879905-umi1 0a45497c6bfbfb22   0
#> 4606209975865882-umi2  0a45497c6bfbfb22   0

# Subset
cg_small <- subset(cg, nodes = CellGraphData(cg, slot = "nodes")[1:100])
cg_small
#> A CellGraph object containing a bipartite graph with 100 nodes and 53 edges
#> Number of markers:  149 
```
