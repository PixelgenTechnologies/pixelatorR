# Find abundant markers by isotype-relative and/or absolute abundance

Identifies markers that are positive in a sufficient fraction of cells.
A cell is positive for a marker if that cell's counts-per-million (CPM)
clears **every** cutoff that is supplied. CPM is computed from counts
plus a pseudocount of 1, so no entry is zero: \\\mathrm{CPM}\_{g,c} =
10^6 (x\_{g,c} + 1) / \sum\_{g'}(x\_{g',c} + 1)\\, where \\x\_{g,c}\\ is
the count for marker \\g\\ in cell \\c\\.

## Usage

``` r
FindAbundantMarkers(
  object,
  isotype_markers,
  isotype_ratio = 1.5,
  abundance_threshold = NULL,
  min_cell_fraction = 0.05,
  group_column = NULL,
  return_stats = FALSE
)
```

## Arguments

- object:

  A `Seurat` object with a `counts` layer. Call
  [`JoinLayers()`](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
  first if the object was merged and still has split layers.

- isotype_markers:

  Character vector of isotype control marker names (for example
  `c("mIgG1", "mIgG2a", "mIgG2b")`).

- isotype_ratio:

  Numeric relative cutoff versus the median isotype CPM, or `NULL` to
  skip this cutoff. Default is `1.5`.

- abundance_threshold:

  Numeric absolute CPM cutoff, or `NULL` to skip this cutoff. Default is
  `NULL`.

- min_cell_fraction:

  Minimum fraction of cells that must be positive for a marker to be
  kept. Default is `0.05`, small enough to keep markers that are
  abundant only in a minority population of a mixed sample.

- group_column:

  Optional metadata column name (typically cell type). If provided, the
  isotype median is computed within each group and the function returns
  a named list with one result per observed group.

- return_stats:

  Logical; if `TRUE`, return a tibble of per-marker statistics instead
  of marker names. The tibble includes `median_cpm`, the median CPM of
  each marker across the cells in the call (within the group when
  `group_column` is set), which shows the CPM range to use when choosing
  `abundance_threshold`. Isotype controls are included.
  `isotype_median_cpm` is the median across all isotype values, not the
  per-marker median. A `kept` column records the filter result.

## Value

If `group_column` is `NULL`: a character vector of kept marker names, or
a tibble if `return_stats` is `TRUE`. If `group_column` is set: a named
list of the same, one element per group.

## Details

- If `isotype_ratio` is set, CPM must be greater than `isotype_ratio`
  times the median CPM of `isotype_markers`.

- If `abundance_threshold` is set, CPM must be greater than
  `abundance_threshold`.

- If both are set, both conditions must hold (AND). Neither argument
  overrides the other. At least one cutoff is required.

A marker is kept if at least `min_cell_fraction` of cells are positive.
Isotype controls are always dropped from the result.

The isotype median is computed from the cells that enter a given call to
the filter. Without `group_column`, that is all cells in `object`. With
`group_column` (for example cell type), the median, positivity, and
`min_cell_fraction` rule are computed independently in each group so a
marker that is absent in one type does not pull the cutoff down for
another. Unused factor levels are omitted.

`isotype_ratio` defaults to `1.5`, matching the relative cutoff used by
the original marker-filter helper. `min_cell_fraction` defaults to
`0.05` so markers that are abundant only in a small population can still
be kept when the function is applied to a whole sample.

## Examples

``` r
library(pixelatorR)

seur <-
  ReadPNA_Seurat(
    minimal_pna_pxl_file(),
    overwrite = TRUE,
    load_proximity_scores = FALSE,
    verbose = FALSE
  )
#> duckdb keeps downloaded extensions and secrets in a temporary directory:
#> ℹ /tmp/Rtmp3BjCtL/duckdb
#> This is removed when the R session ends.
#> • Extensions are re-downloaded each session.
#> • Secrets are lost.
#> ℹ Run duckdb(shared_home = TRUE) (or create ~/.duckdb) to keep them (suitable for most users).
#> ℹ Run duckdb(shared_home = FALSE) to accept the temporary directory (and silence this message).
#> ℹ See ?duckdb_storage for details and alternatives.

# Keep markers above 1.5 times the isotype median CPM in at least 5% of cells
kept <- FindAbundantMarkers(
  object = seur,
  isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b")
)

# Require both the isotype ratio and an absolute CPM floor
kept_and <- FindAbundantMarkers(
  object = seur,
  isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
  isotype_ratio = 1.5,
  abundance_threshold = 2000,
  min_cell_fraction = 0.05
)

# Per-group filtering (isotype median computed within each group)
seur$sample <- c("S1", "S1", "S2", "S2", "S2")
kept_by_sample <- FindAbundantMarkers(
  object = seur,
  isotype_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
  group_column = "sample"
)
```
