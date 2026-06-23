# Load CellGraphs

Loads a list of [`CellGraph`](CellGraph-class.md) objects.

## Usage

``` r
LoadCellGraphs(object, ...)

# S3 method for class 'FileSystemDataset'
LoadCellGraphs(
  object,
  cells,
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  data_type = c("PNA", "MPX"),
  chunk_size = 10,
  verbose = TRUE,
  ...
)

# S3 method for class 'MPXAssay'
LoadCellGraphs(
  object,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'CellGraphAssay'
LoadCellGraphs(
  object,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'CellGraphAssay5'
LoadCellGraphs(
  object,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'PNAAssay'
LoadCellGraphs(
  object,
  cells = colnames(object),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'PNAAssay5'
LoadCellGraphs(
  object,
  cells = colnames(object),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'Seurat'
LoadCellGraphs(
  object,
  assay = NULL,
  cells = colnames(object),
  load_as = c("bipartite", "Anode", "linegraph"),
  add_marker_counts = TRUE,
  add_layouts = FALSE,
  force = FALSE,
  chunk_size = 10,
  cl = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `Seurat` object or an [`CellGraphAssay`](CellGraphAssay-class.md)
  object

- ...:

  Parameters passed to other methods

- cells:

  A character vector of cell names to load CellGraphs for

- load_as:

  Choose how the cell graph should be represented. This option has no
  effect when `data_type = "PNA"` (see details below)

- add_marker_counts:

  Should marker counts be added to the CellGraph objects? Note that this
  option is ignored for the `data.frame` method when
  `data_type = "PNA"`.

- data_type:

  One of "PNA" or "MPX"

- chunk_size:

  Length of chunks used to load CellGraphs from edge list.

- verbose:

  Print messages

- add_layouts:

  Load layouts from the PXL file if available.

- force:

  Force load graph(s) if they are already loaded

- cl:

  A cluster object created by makeCluster, or an integer to indicate
  number of child-processes (integer values are ignored on Windows) for
  parallel evaluations. See Details on performance in the documentation
  for `pbapply`. The default is NULL, which means that no
  parallelization is used.

- assay:

  Assay name

## Value

An object with a list of [`CellGraph`](CellGraph-class.md) objects

## Details

`LoadCellGraphs` will only work if the PXL file path(s) stored in the
object are valid. You can check the PXL file paths stored in the object
with the [`FSMap`](FSMap.md) method. If the PXL file(s) are invalid, an
error will be thrown. Visit [`RestorePaths`](RestorePaths.md) for more
information on how to update the PXL file paths.

## Graph representations

Note that PNA component graphs are always loaded as bipartite graphs.

For MPX data, one of the following graph representations can be used by
specifying the `load_as` parameter:

In the bipartite graph, edges can only go from a upia to a upib. The
bipartite graph is first collapsed from a multigraph to a simple graph
by aggregating marker counts for each upia/upib combination. For
visualization of marker counts on the graph, it is often convenient to
project values on the nodes; however, the marker counts are not
available in the nodes in the bipartite graph. To circumvent this issue,
node counts are calculated by aggregating its edge counts. This means
that the total marker count will be inflated.

In the A node projected-graph, pairs of upias that share a upib are
connected with an edge. The (upia) node counts are calculated by
aggregating counts across all edges associated with the A nodes. The
number of nodes in the A node projected-graph is substantially lower
than the bipartite graph.

Starting with a bipartite graph, a node is placed on each edge. Then, an
edge is drawn between each pair of adjacent nodes. The number of nodes
in the linegraph is equivalent to the number of edges in the bipartite
graph. However, the number of edges is substantially larger. One major
advantage with the linegraph is that the node counts represent the
actual raw data, which is not the case for the bipartite graph and the A
node projected-graph. On the down side, linegraphs are much larger and
tends to slow down layout computations and visualizations.

## Examples

``` r
# Read from PNAAssay
pxl_file <- minimal_pna_pxl_file()
seur_obj <- ReadPNA_Seurat(pxl_file)
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
pna_assay <- LoadCellGraphs(seur_obj[["PNA"]], cells = "0a45497c6bfbfb22")
#> ℹ Fetching edgelists for 1 cells 
#> → Creating <CellGraph> objects
#> → Fetching marker counts
#> → Adding marker counts to <CellGraph> object(s)
#> ✔ Successfully loaded 1 <CellGraph> object(s).
CellGraphs(pna_assay)[["0a45497c6bfbfb22"]]
#> A CellGraph object containing a bipartite graph with 43543 nodes and 97014 edges
#> Number of markers:  149 
```
