# Join component-level metadata onto a proximity or edgelist table

S3 generic used to join a local `data.frame` of component-level metadata
or marker counts onto an eager or lazy proximity/edgelist table.

## Usage

``` r
join_proximity_metadata(proximity_scores, mdata, by = "component", ...)

# Default S3 method
join_proximity_metadata(proximity_scores, mdata, by = "component", ...)

# S3 method for class 'tbl_lazy'
join_proximity_metadata(
  proximity_scores,
  mdata,
  by = "component",
  name = NULL,
  reuse = FALSE,
  ...
)
```

## Arguments

- proximity_scores:

  A `data.frame` or `tbl_lazy` table.

- mdata:

  A local `data.frame` containing `component` and metadata columns to
  join.

- by:

  Join specification passed to
  [`left_join`](https://dplyr.tidyverse.org/reference/mutate-joins.html).
  Defaults to `"component"`.

- ...:

  Additional arguments passed to methods.

- name:

  Optional name for the temporary table created in the database. When
  `NULL` (default), a unique temporary table name is generated.

- reuse:

  Logical indicating whether to reuse an existing table named `name` on
  the database connection instead of re-uploading `mdata`.

## Value

`proximity_scores` with columns from `mdata` joined.
