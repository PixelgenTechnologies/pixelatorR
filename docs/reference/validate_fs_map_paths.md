# Validate PXL file paths in an `fs_map` table

S3 generic used to validate the `pxl_file` column of an `fs_map` table.
The default method checks that each path exists on the local file
system.

## Usage

``` r
validate_fs_map_paths(fs_map, call = caller_env(), ...)

# Default S3 method
validate_fs_map_paths(fs_map, call = caller_env(), ...)
```

## Arguments

- fs_map:

  A `tbl_df` with columns `id_map`, `sample`, and `pxl_file`.

- call:

  The calling environment for error reporting.

- ...:

  Additional arguments passed to methods.

## Value

`NULL`, invisibly.
