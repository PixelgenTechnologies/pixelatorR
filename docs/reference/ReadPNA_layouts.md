# Load layouts from a PNA PXL file

Load layouts from a PNA PXL file

## Usage

``` r
ReadPNA_layouts(
  pxl_file,
  cells = NULL,
  add_marker_counts = FALSE,
  verbose = FALSE
)
```

## Arguments

- pxl_file:

  Path to a PXL file containing PNA data

- cells:

  A character vector with component IDs. If NULL, all components are
  loaded.

- add_marker_counts:

  Logical specifying if marker counts should be added to the layout
  table(s).

- verbose:

  Logical specifying if verbose output should be printed

## Value

A list with one `tbl_df` object for each component. Each object contains
the node IDs and their x, y, and z coordinates.
