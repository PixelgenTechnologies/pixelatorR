# Lazy load a table from one or more PXL sources

S3 generic used by [`ProximityScores`](ProximityScores.md) and
[`Edgelists`](Edgelists.md) to lazily load tables from the sources
referenced in `fs_map`.

## Usage

``` r
lazy_load_table(
  fs_map,
  table_name = "proximity",
  calc_log2ratio = TRUE,
  proteins_keep = NULL,
  union = TRUE,
  call = caller_env(),
  ...
)

# Default S3 method
lazy_load_table(
  fs_map,
  table_name = "proximity",
  calc_log2ratio = TRUE,
  proteins_keep = NULL,
  union = TRUE,
  call = caller_env(),
  ...
)
```

## Arguments

- fs_map:

  A `tbl_df` with information about the PXL sources and component IDs.

- table_name:

  The name of the table to load, e.g. `"proximity"` or `"edgelist"`.

- calc_log2ratio:

  A logical indicating whether to calculate the log2 ratio proximity
  score. Only used for the `"proximity"` table.

- proteins_keep:

  A character vector of protein names (`marker_1`/`marker_2`) to keep in
  the proximity score table.

- union:

  A logical indicating whether to union the results from multiple
  sources into a single table. If `FALSE`, a named list of lazy tables
  (one per source) is returned.

- call:

  The calling environment for error reporting.

- ...:

  Additional arguments passed to methods.

## Value

A `tbl_lazy` table (when `union = TRUE`) or a named list of `tbl_lazy`
tables (when `union = FALSE`).
