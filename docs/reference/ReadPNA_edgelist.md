# Load the edge list from a PNA PXL file

Note that the umi1 and umi2 sequences are encoded as `integer64` which
is not natively supported by R. The `bit64` package is required to read
these IDs and conversion to `integer` or `numeric` should be avoided. It
is safer to convert these to character vectors if needed for downstream
processing.

## Usage

``` r
ReadPNA_edgelist(
  pxl_file,
  cells = NULL,
  umi_data_type = c("int64", "string", "suffixed_string"),
  lazy = TRUE
)
```

## Arguments

- pxl_file:

  Path to a PXL file with a PNA data set

- cells:

  A character vector with component IDs. If NULL, all components are
  loaded.

- umi_data_type:

  One of "int64", "string" or "suffixed_string". Default is "int64".

  - "int64": The UMIs are encoded as int64

  - "string": The UMIs are encoded as character

  - "suffixed_string": The UMIs are encoded as character with a suffix
    '-umi1' or '-umi2' added

- lazy:

  A logical specifying whether to load the data lazily. If `TRUE`, a
  `tbl_lazy` object is returned.

## Value

A `tbl_df` or, if `lazy = TRUE`, a `tbl_lazy` object with the edge list:

- marker_1: The first protein

- marker_2: The second protein

- umi1: A unique ID of the first RCA product

- umi2: A unique ID of the second RCA product

- read_count: The number of reads supporting the edge

- uei_count: The number of unique event identifiers (UEIs) supporting
  the edge. This column is optional and is only returned if it is
  present in the edgelist table.

- component: The component name
