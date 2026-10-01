# Create a NodeDimReduc object

Constructs a [`NodeDimReduc`](NodeDimReduc-class.md) for storage in the
`reductions` slot of a [`CellGraph`](CellGraph-class.md). Analogous to
[`CreateDimReducObject`](https://satijalab.github.io/seurat-object/reference/CreateDimReducObject.html).

## Usage

``` r
CreateNodeDimReducObject(
  embeddings,
  loadings = NULL,
  stdev = numeric(),
  key = "DR_",
  method = character(),
  misc = list()
)
```

## Arguments

- embeddings:

  A numeric matrix of node embeddings. Row names must be node names. If
  `colnames` are missing they are set to
  `paste0(key, seq_len(ncol(embeddings)))`.

- loadings:

  An optional numeric matrix of feature loadings. Columns are aligned to
  the embedding dimensions.

- stdev:

  An optional numeric vector of standard deviations, one per dimension

- key:

  Prefix for embedding column names. A trailing underscore is added if
  missing.

- method:

  Optional name of the reduction method

- misc:

  A list of additional metadata to store with the reduction

## Value

A [`NodeDimReduc`](NodeDimReduc-class.md) object
