# Differential analysis (polarity)

Runs differential analysis on polarity scores generated with the
`Pixelator` data processing pipeline.

## Usage

``` r
RunDPA(object, ...)

# S3 method for class 'data.frame'
RunDPA(
  object,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  polarity_metric = "morans_z",
  min_n_obs = 0,
  cl = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
)

# S3 method for class 'Seurat'
RunDPA(
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  polarity_metric = "morans_z",
  min_n_obs = 0,
  cl = NULL,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object containing polarity scores

- ...:

  Not yet implemented

- contrast_column:

  The name of the column where the group labels are stored. This column
  must include `target` and `reference`.

- reference:

  The name of the reference group

- targets:

  The names of the target groups. These groups will be compared to the
  reference group. If the value is set to `NULL` (default), all groups
  available in `contrast_column` except `reference` will be compared to
  the `reference` group.

- group_vars:

  An optional character vector with column names to group the tests by.

- polarity_metric:

  The polarity metric to use. Any numeric data column in the polarity
  score table can be selected. The default is "morans_z".

- min_n_obs:

  Minimum number of observations allowed in a group. Target groups with
  less observations than `min_n_obs` will be skipped.

- cl:

  A cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html), or an
  integer to indicate number of child-processes (integer values are
  ignored on Windows) for parallel evaluations. See Details on
  performance in the documentation for `pbapply`. The default is NULL,
  which means that no parallelization is used. Note that warnings are
  not caught when using parallel processing.

- alternative:

  One of 'two.sided', 'less' or 'greater' (see
  [`?wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html) for
  details)

- conf_int:

  Should confidence intervals be computed? (see
  [`?wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html) for
  details)

- p_adjust_method:

  One of "bonferroni", "holm", "hochberg", "hommel", "BH", "BY" or
  "fdr". (see [`?p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for
  details)

- verbose:

  Print messages

- assay:

  Name of assay to use

## Value

A `tbl_df` object with test results

## Details

If you are working with a `Seurat` object created with pixelatorR that
contains a [`CellGraphAssay`](CellGraphAssay-class.md), the polarity
scores are accessed directly from the
[`CellGraphAssay`](CellGraphAssay-class.md) (see
[`PolarizationScores`](PolarizationScores.md)).

The input object should contain a `contrast_column` (character vector or
factor) that includes information about the groups to compare. A typical
example is a column with sample labels, for instance: "control",
"stimulated1", "stimulated2". If the input object is a `Seurat` object,
the `contrast_column` should be available in the `meta.data` slot. For
those familiar with `FindMarkers` from Seurat, `contrast_column` is
equivalent to the `group.by` parameter.

The `targets` parameter specifies a character vector with the names of
the groups to compare `reference`. `targets` can be a single group name
or a vector of group names while `reference` can only refer to a single
group. Both `targets` and `reference` should be present in the
`contrast_column`. These parameters are similar to the `ident.1` and
`ident.2` parameters in `FindMarkers`.

## Additional groups

The test is always computed between `targets` and `reference`, but it is
possible to add additional grouping variables with `group_vars`. If
`group_vars` is used, each comparison is split into groups defined by
the `group_vars`. For instance, if we have annotated cells into cell
type populations and saved these annotations in a `meta.data` column
called "cell_type", we can pass "cell_type" to `group_vars="cell_type"`
to split tests across each cell type.

## Types of comparisons

Consider a scenario where we have a Seurat object (`seurat_object`) with
PNA data. `seurat_object` contains a `meta.data` column called
"sampleID" that holds information about what samples the PNA components
originated from. This column could have three sample IDs: "control",
"stimulated1" and "stimulated2". In addition, we have a column called
"cell_type" that holds information about the cell type identity of each
PNA component.

1.  If we want to compare the "stimulated1" group to the "control"
    group:

        dpa_markers <- RunDPA(object = seurat_object,
                              contrast_column = "sampleID",
                              reference = "control",
                              targets = "stimulated1")

2.  If we want to compare the "stimulated1" and "stimulated2" groups to
    the "control" group:

        dpa_markers <- RunDPA(object = seurat_object,
                              contrast_column = "sampleID",
                              reference = "control",
                              targets = c("stimulated1", "stimulated2"))

3.  If we want to compare the "stimulated1" and "stimulated2" groups to
    the "control" group, and split the tests by cell type:

        dpa_markers <- RunDPA(object = seurat_object,
                             contrast_column = "sampleID",
                             reference = "control",
                             targets = c("stimulated1", "stimulated2"),
                             group_vars = "cell_type")

## Error handling

If the test fails for a certain comparison, a warning is raised and no
results will be returned for that comparison. This can happen if one of
the two groups being compared has too few observations. Note that when
using parallel processing, these warnings are muted to reduce the
overhead of communication between processes.

## See also

Other DA-methods:
[`DifferentialProximityAnalysis()`](DifferentialProximityAnalysis.md),
[`RunDAA()`](RunDAA.md), [`RunDCA()`](RunDCA.md)

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_mpx_pxl_file()

# Load polarization scores
polarization_table1 <- polarization_table2 <- ReadMPX_polarization(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de0b91932
polarization_table1$sample <- "Sample1"
polarization_table2$sample <- "Sample2"
polarization_table_merged <- bind_rows(polarization_table1, polarization_table2)

# Run DPA using table as input
dpa_markers <- RunDPA(polarization_table_merged,
  contrast_column = "sample",
  targets = "Sample1", reference = "Sample2"
)
#> Warning: Got the following message when running wilcox.test test for marker 'ACTB': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'B2M': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD102': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD11a': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD11b': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD11c': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD127': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD137': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD14': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD150': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD152': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD154': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD158': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD16': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD161': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD162': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD163': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD18': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD19': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD197': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD1d': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD2': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD20': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD200': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD22': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD229': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD244': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD25': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD26': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD268': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD27': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD274': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD278': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD279': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD29': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD314': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD32': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD328': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD33': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD337': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD35': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD36': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD37': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD38': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD3E': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD4': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD40': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD41': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD43': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD44': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45RA': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45RB': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD47': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD48': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD49D': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD5': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD50': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD52': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD53': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD54': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD55': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD59': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD62P': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD64': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD69': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD7': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD71': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD72': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD8': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD82': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD84': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD86': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD9': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'HLA-ABC': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'HLA-DR': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'TCRb': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'mIgG1': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'mIgG2a': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'mIgG2b': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
dpa_markers
#> # A tibble: 80 × 14
#>      estimate data_type target  reference    n1    n2 statistic     p p_adj
#>         <dbl> <chr>     <chr>   <chr>     <int> <int>     <dbl> <dbl> <dbl>
#>  1  0         morans_z  Sample1 Sample2       4     4       8       1     1
#>  2  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  3  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  4 -0.0000610 morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  5  0         morans_z  Sample1 Sample2       4     4       8       1     1
#>  6  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  7  0.0000256 morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  8  0         morans_z  Sample1 Sample2       3     3       4.5     1     1
#>  9  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#> 10  0.0000178 morans_z  Sample1 Sample2       5     5      12.5     1     1
#> # ℹ 70 more rows
#> # ℹ 5 more variables: conf.low <dbl>, conf.high <dbl>, method <chr>,
#> #   alternative <chr>, marker <chr>

# Seurat objects
seur1 <- seur2 <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de019676db2
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/RtmpInAYcT/dir7de0162c1f1f
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/RtmpInAYcT/file7de06409662f.h5ad
seur1$sample <- "Sample1"
seur2$sample <- "Sample2"
seur_merged <- merge(seur1, seur2, add.cell.ids = c("A", "B"))

# Run DPA
dpa_markers <- RunDPA(seur_merged,
  contrast_column = "sample",
  targets = "Sample1", reference = "Sample2"
)
#> Warning: Got the following message when running wilcox.test test for marker 'ACTB': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'B2M': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD102': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD11a': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD11b': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD11c': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD127': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD137': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD14': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD150': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD152': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD154': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD158': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD16': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD161': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD162': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD163': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD18': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD19': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD197': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD1d': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD2': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD20': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD200': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD22': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD229': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD244': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD25': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD26': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD268': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD27': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD274': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD278': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD279': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD29': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD314': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD32': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD328': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD33': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD337': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD35': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD36': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD37': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD38': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD3E': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD4': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD40': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD41': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD43': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD44': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45RA': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45RB': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD47': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD48': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD49D': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD5': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD50': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD52': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD53': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD54': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD55': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD59': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD62P': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD64': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD69': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD7': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD71': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD72': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD8': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD82': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD84': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD86': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD9': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'HLA-ABC': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'HLA-DR': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'TCRb': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'mIgG1': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'mIgG2a': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'mIgG2b': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
dpa_markers
#> # A tibble: 80 × 14
#>      estimate data_type target  reference    n1    n2 statistic     p p_adj
#>         <dbl> <chr>     <chr>   <chr>     <int> <int>     <dbl> <dbl> <dbl>
#>  1  0         morans_z  Sample1 Sample2       4     4       8       1     1
#>  2  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  3  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  4 -0.0000610 morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  5  0         morans_z  Sample1 Sample2       4     4       8       1     1
#>  6  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  7  0.0000256 morans_z  Sample1 Sample2       5     5      12.5     1     1
#>  8  0         morans_z  Sample1 Sample2       3     3       4.5     1     1
#>  9  0         morans_z  Sample1 Sample2       5     5      12.5     1     1
#> 10  0.0000178 morans_z  Sample1 Sample2       5     5      12.5     1     1
#> # ℹ 70 more rows
#> # ℹ 5 more variables: conf.low <dbl>, conf.high <dbl>, method <chr>,
#> #   alternative <chr>, marker <chr>
```
