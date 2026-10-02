# The NodeDimReduc class

The `NodeDimReduc` class stores a dimensionality reduction computed on
the nodes of a [`CellGraph`](CellGraph-class.md). It follows the same
design as
[`DimReduc`](https://satijalab.github.io/seurat-object/reference/DimReduc-class.html):
a required embeddings matrix, optional feature loadings, standard
deviations, a dimension key, and a miscellaneous list for
method-specific extras.

## Slots

- `embeddings`:

  A numeric `matrix` of node embeddings (nodes x dimensions). Row names
  are node names when the object is created. When the reduction is
  stored on a [`CellGraph`](CellGraph-class.md), embedding rows follow
  the graph `nodes` map and row names are dropped.

- `loadings`:

  An optional numeric `matrix` of feature loadings (features x
  dimensions)

- `stdev`:

  A numeric vector of standard deviations (or eigenvalues) for each
  dimension

- `key`:

  A character scalar used as the column-name prefix, ending with `_`
  (for example `"PC_"`)

- `method`:

  A character scalar naming the reduction method (for example `"pca"` or
  `"umap"`)

- `misc`:

  A named list of unstructured additional data
