# Load layouts from an PXL file containing MPX data

Layouts can be pre-computed with the Pixelator data processing pipeline
and are stored in a hive-styled partitioning in the PXL file. This
function reads the layouts from the PXL file and returns them as a list.
Use [`inspect_pxl_file`](inspect_pxl_file.md) to check the contents of a
PXL file.

## Usage

``` r
ReadMPX_layouts(
  filename,
  cells = NULL,
  graph_projection = c("bipartite", "Anode", "linegraph", "full"),
  verbose = TRUE
)
```

## Arguments

- filename:

  Path to a PXL file

- cells:

  A character vector with component IDs. If NULL, all components are
  loaded.

- graph_projection:

  The graph projection to load. Default is 'bipartite'. If multiple
  projections are present in the file, only the selected one is loaded.

- verbose:

  Print messages

## Value

A list of lists with the layouts. At the top level, the list is split by
layout. At the second level, the list is split by component. The
components are sorted in the order they appear in the PXL file.
