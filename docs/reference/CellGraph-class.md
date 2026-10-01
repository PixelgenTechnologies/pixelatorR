# The CellGraph class

The CellGraph class is designed to hold information needed for working
with PNA single-cell graphs.

## Details

A `CellGraph` contains counts and a graph, and optionally layouts,
layers, metadata, or reductions.

Node-level variable names must be unique across the graph node table,
`meta.data`, and reduction embeddings, and must not overlap count or
layer features. Count and layer matrices are the exception: they may
share feature names because callers select a layer explicitly.

## Slots

- `cellgraph`:

  A `tbl_graph` object corresponding to a cell graph

- `nodes`:

  Character vector of node IDs in graph order. This is the map used to
  align counts, layouts, layers, metadata, and reductions. Those tables
  are stored in this order without copying the IDs as row names.

- `counts`:

  A `matrix`-like object with marker counts (nodes x markers). Rows
  follow `nodes`. The counts matrix can be extracted as the `"counts"`
  layer via
  [`Layers`](https://satijalab.github.io/seurat-object/reference/Layers.html)
  /
  [`LayerData`](https://satijalab.github.io/seurat-object/reference/Layers.html).

- `layout`:

  A named `list` of `data.frame` objects with coordinates for cell
  layouts. Rows follow `nodes`. A `name` column or explicit row names
  are accepted on input and used only to reorder; MPX bipartite layouts
  may omit `-A`/`-B` suffixes. Stored layouts keep coordinate columns
  (typically `x`, `y`, `z`). Layouts without node IDs still work if the
  number of rows matches the graph.

- `layers`:

  A named `list` of additional numeric node matrices (nodes x features),
  analogous to layers on a Seurat
  [`Assay5`](https://satijalab.github.io/seurat-object/reference/Assay5-class.html).
  A layer can be extracted via
  [`Layers`](https://satijalab.github.io/seurat-object/reference/Layers.html)
  /
  [`LayerData`](https://satijalab.github.io/seurat-object/reference/Layers.html).

- `meta.data`:

  A `data.frame` of node-level metadata (one row per node). Rows follow
  `nodes`. Columns may have mixed types.

- `reductions`:

  A named `list` of [`NodeDimReduc`](NodeDimReduc-class.md) objects
