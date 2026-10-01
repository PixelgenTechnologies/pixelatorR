# Compute local proximity scores

Computes local proximity scores (LPS) for nodes in PNA cell graphs using
[`local_proximity`](local_proximity.md) and stores the result on each
[`CellGraph`](CellGraph-class.md).

Computes local proximity scores for nodes in PNA cell graphs using
[`local_proximity`](local_proximity.md) and stores the output on each
[`CellGraph`](CellGraph-class.md).

## Usage

``` r
ComputeLPS(object, ...)

# S3 method for class 'CellGraph'
ComputeLPS(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  A_k = NULL,
  seed = 123,
  name = "lps",
  ...
)

# S3 method for class 'CellGraphList'
ComputeLPS(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
)

# S3 method for class 'PNAAssay'
ComputeLPS(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
)

# S3 method for class 'PNAAssay5'
ComputeLPS(
  object,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
)

# S3 method for class 'Seurat'
ComputeLPS(
  object,
  assay = NULL,
  markers = NULL,
  method = c("analytical", "permutation"),
  mode = c("self-clustering", "all", "any"),
  iterations = 50L,
  k = 3L,
  seed = 123,
  name = "lps",
  verbose = TRUE,
  cl = NULL,
  ...
)
```

## Arguments

- object:

  An object

- ...:

  Additional parameters passed to other methods

- markers:

  A character vector specifying the markers to use. If `NULL`, all
  markers in the count matrix of each `CellGraph` are used. Methods that
  iterate over multiple `CellGraph` objects keep the intersection with
  available markers. If a graph has no counts or none of the requested
  markers, a warning is emitted and that graph is left unmodified. The
  single-graph method still errors when counts are missing.

- method:

  A character string specifying the method to use for computing the
  local proximity score. Options are `"analytical"` or `"permutation"`.

- mode:

  A character string specifying the mode of computation. See
  [`local_proximity`](local_proximity.md) for details. Default is
  `"self-clustering"`.

- iterations:

  An integer specifying the number of iterations to run when
  `method = "permutation"`.

- k:

  An integer specifying the neighborhood size to consider.

- A_k:

  An optional pre-computed expanded adjacency matrix. Only used by the
  `CellGraph` method.

- seed:

  An integer for random seed setting.

- name:

  Name of the layer (when a matrix is returned) or metadata column (when
  a vector is returned) used to store scores. Default is `"lps"`.

- verbose:

  Print messages

- cl:

  An integer to indicate number of child-processes (integer values are
  ignored on Windows) for parallel evaluations. See Details on
  performance in the documentation for `pbapply`. The default is `NULL`,
  which means that no parallelization is used.

- assay:

  Name of PNAAssay containing the cell graphs to compute local proximity
  scores

## Value

An object with local proximity scores stored on each
[`CellGraph`](CellGraph-class.md). Matrix results are stored as a layer;
vector results are stored in node `meta.data`.

## Details

The default `mode` is `"self-clustering"`, which returns a
node-by-marker matrix that is stored as a layer. Other modes return a
single score per node, which is stored in `meta.data`.

## Examples

``` r
library(pixelatorR)

se <- ReadPNA_Seurat(minimal_pna_pxl_file(), verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpmC3mql/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
se <- LoadCellGraphs(se, cells = colnames(se)[1], verbose = FALSE)
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/RtmpmC3mql/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.
cg <- CellGraphs(se)[[1]]

# Matrix result is stored as a layer
cg <- ComputeLPS(cg, markers = "B2M")
#> Warning: 'as(<dgCMatrix>, "ngCMatrix")' is deprecated.
#> Use 'as(., "nMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
SeuratObject::Layers(cg)
#> [1] "counts" "lps"   

# Vector result is stored in node metadata
cg <- ComputeLPS(cg, markers = "B2M", mode = "all", name = "lps_b2m")
head(CellGraphData(cg, slot = "meta.data"))
#>      lps_b2m
#> 1  0.3510601
#> 2 -1.3058178
#> 3  0.8001271
#> 4  1.2880565
#> 5  1.5179853
#> 6  1.4691879

# Compute LPS for loaded cell graphs in a PNAAssay
pna_assay <- ComputeLPS(se[["PNA"]], markers = "B2M")
#> ℹ Computing local proximity scores for 1 graph

# Seurat method (only loaded cell graphs are processed)
seur <- ComputeLPS(se, markers = "B2M")
#> ℹ Computing local proximity scores for 1 graph
```
