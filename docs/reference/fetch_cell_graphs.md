# Fetch cell graphs from PXL sources

S3 generic used by [`LoadCellGraphs`](LoadCellGraphs.md) on
[`PNAAssay`](PNAAssay-class.md) and [`PNAAssay5`](PNAAssay5-class.md)
objects to extract [`CellGraph`](CellGraph-class.md) objects for a set
of cells from the sources referenced in `fs_map`.

## Usage

``` r
fetch_cell_graphs(
  fs_map,
  cells_to_load,
  add_layouts = FALSE,
  verbose = TRUE,
  ...
)

# Default S3 method
fetch_cell_graphs(
  fs_map,
  cells_to_load,
  add_layouts = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- fs_map:

  A `tbl_df` with columns `id_map`, `sample`, and `pxl_file`.

- cells_to_load:

  A character vector of current cell IDs to load.

- add_layouts:

  Logical indicating whether to load precomputed layouts.

- verbose:

  Print messages.

- ...:

  Additional arguments passed to methods.

## Value

A named list of [`CellGraph`](CellGraph-class.md) objects keyed by
`cells_to_load`.
