# Filter an edgelist by downsampling read counts

Filter an edgelist (for a subset of components) in a PXL file by
downsampling read counts using predefined fractions. The probability of
a read staying in the edgelist is \\1 - (1 - frac) ^ read_count\\.

## Usage

``` r
downsample_to_parquet(pxl_file, outdir, components, fracs = seq(0.1, 0.9, 0.1))
```

## Arguments

- pxl_file:

  Path to the PXL file.

- outdir:

  Directory to write the output parquet files.

- components:

  A character vector of component names to keep.

- fracs:

  A numeric vector of fractions to downsample by.

## Value

A tibble with paths to the written parquet files and the fractions.

## Details

The filtered edgelist(s) are written to parquet files in `outdir`.
