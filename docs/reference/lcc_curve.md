# Compute LCC curve for downsampled edgelists

**\[experimental\]**

Computes the sizes of the largest connected components (LCC) for
edgelists downsampled at multiple fractions. This is useful for
analyzing graph stability and determining the critical sequencing depth.

## Usage

``` r
lcc_curve(
  pxl_file,
  components,
  fracs = seq(0.1, 1, by = 0.1),
  outdir = NULL,
  mc_cores = 1,
  verbose = TRUE
)
```

## Arguments

- pxl_file:

  Path to the input PXL file containing the edgelist.

- components:

  A character vector of component names to include. All components must
  be present in the PXL file.

- fracs:

  A numeric vector of downsampling fractions. Values must be strictly
  between 0 and 1. Default is `seq(0.1, 1, by = 0.1)`.

- outdir:

  Optional output directory for downsampled parquet files. If `NULL`
  (default), files are saved to a temporary directory and deleted after
  computation.

- mc_cores:

  Number of cores for parallel processing. Default is 1 (sequential).
  Values \> 1 use
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html).

- verbose:

  Logical; if `TRUE` (default), print progress messages.

## Value

A tibble with columns:

- component:

  The component (cell) identifier.

- frac:

  The downsampling fraction.

- n_nodes:

  The number of nodes in the largest connected component.

- read_count:

  The estimated number of reads supporting the LCC at the given
  fraction.

## Details

### Graph stability analysis

The LCC (largest connected component) measures graph connectivity. When
downsampling sequencing reads, edges are removed probabilistically,
which eventually causes the graph to fragment into smaller disconnected
components.

By computing LCC sizes at multiple downsampling fractions, you can
observe how graph connectivity degrades with reduced sequencing depth.
Typically, there is a sharp transition point where the LCC collapses,
indicating the minimum sequencing depth required to maintain graph
structure.

### Workflow

This function is a high-level wrapper that:

1.  Calls [`downsample_to_parquet`](downsample_to_parquet.md) to create
    downsampled edgelists

2.  Calls [`lcc_sizes`](lcc_sizes.md) to compute LCC sizes for each
    fraction

3.  Optionally cleans up temporary files

### Performance

Computations are scalable thanks to efficient SQL queries. LCCs are
computed using the `duckpgq` extension for DuckDB, enabling fast graph
operations without loading the full edgelist into memory.

For large datasets, it is recommended to compute LCC sizes on a subset
of components to minimize memory usage and computation time.

## See also

[`downsample_to_parquet`](downsample_to_parquet.md) for the downsampling
step. [`lcc_sizes`](lcc_sizes.md) for computing LCC sizes from parquet
files. [`approximate_saturation_curve`](approximate_saturation_curve.md)
for saturation analysis without actual downsampling.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
library(dplyr)
pxl_file <- minimal_pna_pxl_file()
components <- ReadPNA_counts(pxl_file) %>% colnames()

# Select fractions based on average read depth
db <- PixelDB$new(pxl_file)
avg_reads <- db$cell_meta()$reads %>% mean()
fracs <- (10^seq(4, log10(avg_reads), length.out = 20))[-20] / avg_reads

lcc_df <- lcc_curve(pxl_file, components, fracs = fracs)

# Plot LCC sizes by fraction
ggplot(lcc_df, aes(frac, n_nodes, color = component)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  guides(color = "none")

# Compute graph stability (LCC / max theoretical nodes)
nodesat <- approximate_node_saturation(db) %>% collect()
db$close()

lcc_df <- lcc_df %>%
  left_join(nodesat, by = "component")

ggplot(lcc_df, aes(frac, n_nodes / nodes, color = component)) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(y = "Graph stability (LCC / total nodes)") +
  theme_bw() +
  guides(color = "none")
} # }
```
