# Compute LCC sizes for downsampled edgelists

**\[experimental\]** This function computes the sizes of the largest
connected components (LCC) for downsampled edgelists. The downsampling
is determined by the `fracs` parameter, which specifies the fractions of
sequencing reads to retain.

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

  A character vector of component names to include in the computation.
  These must be present in the PXL file.

- fracs:

  A numeric vector specifying the fractions of sequencing reads to
  retain.

- outdir:

  Optional output directory where the downsampled parquet files will be
  saved. If `NULL`, the files will be saved in a temporary directory
  which is deleted post processing.

- mc_cores:

  Number of cores to use for parallel processing. Default is 1
  indicating sequential processing.

- verbose:

  Logical indicating whether to print progress messages.

## Value

A tibble with columns `component` (the original component ID), `frac`
(the fraction) and `n_nodes` (the size of the LCC for that fraction and
component).

## LCC

The LCC (largest connected component) can be used to analyze graph
stability. When downsampling sequencing reads, edges are removed, which
eventually leads to a collapse of the graph into many smaller
components. By downsampling the graph at multiple fractions, we can
observe how the size of the largest connected component changes.
Typically, there is a sharp drop in the size of the LCC at a certain
fraction, indicating that the graph no longer has a large connected
component and becomes completely fragmented. Note that the LCC is
computed for each fraction and component (cell) separately.

## NA

The LCC (largest connected component) can be used to analyze graph
stability. When downsampling sequencing reads, edges are removed, which
eventually leads to a collapse of the graph into many smaller
components. By downsampling the graph at multiple fractions, we can
observe how the size of the largest connected component changes.
Typically, there is a sharp drop in the size of the LCC at a certain
fraction, indicating that the graph no longer has a large connected
component and becomes completely fragmented. Note that the LCC is
computed for each fraction and component (cell) separately.

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)
library(dplyr)
pxl_file <- minimal_pna_pxl_file()
components <- ReadPNA_counts(pxl_file) %>% colnames()

# Select fractions
db <- PixelDB$new(pxl_file)
avg_reds <- db$cell_meta()$reads %>% mean()
fracs <- (10^seq(4, log10(avg_reds), length.out = 20))[-20] / avg_reds

lcc_df <- lcc_curve(
  minimal_pna_pxl_file(),
  components,
  fracs = fracs
)

ggplot(lcc_df, aes(frac, n_nodes, color = component)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  guides(color = "none")

# We can calculate graph stability by dividing
# the LCC by the maximum theoretical number of nodes
nodesat <- approximate_node_saturation(db) %>%
  collect()
db$close()
lcc_df <- lcc_df %>%
  left_join(nodesat, by = "component")

ggplot(lcc_df, aes(frac, n_nodes / theoretical_max_nodes, color = component)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  guides(color = "none") +
  scale_y_continuous(limits = c(0, 1))
} # }
```
