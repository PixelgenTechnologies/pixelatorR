# Differential analysis (colocalization)

Runs differential analysis on colocalization scores generated with the
`Pixelator` data processing pipeline.

## Usage

``` r
RunDCA(object, ...)

# S3 method for class 'data.frame'
RunDCA(
  object,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  coloc_metric = "pearson_z",
  min_n_obs = 0,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  cl = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'Seurat'
RunDCA(
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  coloc_metric = "pearson_z",
  min_n_obs = 0,
  alternative = c("two.sided", "less", "greater"),
  conf_int = TRUE,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  cl = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object containing colocalization scores

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

- coloc_metric:

  The colocalization metric to use. Any numeric data column in the
  colocalization score table can be selected. The default is
  "pearson_z".

- min_n_obs:

  Minimum number of observations allowed in a group. Target groups with
  less observations than `min_n_obs` will be skipped.

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

- cl:

  A cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html), or an
  integer to indicate number of child-processes (integer values are
  ignored on Windows) for parallel evaluations. See Details on
  performance in the documentation for `pbapply`. The default is NULL,
  which means that no parallelization is used. Note that warnings are
  not caught when using parallel processing.

- verbose:

  Print messages

- assay:

  Name of assay to use

## Value

A `tbl_df` object with test results

## Details

If you are working with a `Seurat` object created with pixelatorR that
contains a [`CellGraphAssay`](CellGraphAssay-class.md), the
colocalization scores are accessed directly from the
[`CellGraphAssay`](CellGraphAssay-class.md) (see
[`ColocalizationScores`](ColocalizationScores.md)).

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

        dca_markers <- RunDCA(object = seurat_object,
                              contrast_column = "sampleID",
                              reference = "control",
                              targets = "stimulated1")

2.  If we want to compare the "stimulated1" and "stimulated2" groups to
    the "control" group:

        dca_markers <- RunDCA(object = seurat_object,
                              contrast_column = "sampleID",
                              reference = "control",
                              targets = c("stimulated1", "stimulated2"))

3.  If we want to compare the "stimulated1" and "stimulated2" groups to
    the "control" group, and split the tests by cell type:

        dca_markers <- RunDCA(object = seurat_object,
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
[`RunDAA()`](RunDAA.md), [`RunDPA()`](RunDPA.md)

## Examples

``` r
library(pixelatorR)
library(dplyr)

pxl_file <- minimal_mpx_pxl_file()
# Seurat objects
seur1 <- seur2 <- ReadMPX_Seurat(pxl_file)
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf4ecf3840
#> ! Failed to remove temporary dir C:/Users/max/AppData/Local/Temp/Rtmpampkmn/dir5bf414b06d87
#> ! Failed to remove temporary file C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf469876177.h5ad
seur1$sample <- "Sample1"
seur2$sample <- "Sample2"
seur_merged <- merge(seur1, seur2, add.cell.ids = c("A", "B"))

# Subset data to run test on a few markers
seur_merged <- subset(seur_merged,
  features = c(
    "CD3", "CD4", "CD8", "CD19",
    "CD20", "CD45RA", "CD45RO"
  )
)

# Run DCA
dca_markers <- RunDCA(seur_merged,
  contrast_column = "sample",
  target = "Sample1", reference = "Sample2"
)
#> Warning: Got the following message when running wilcox.test test for marker 'CD19/CD20': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD19/CD4': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD19/CD45RA': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD19/CD8': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD20/CD4': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD20/CD45RA': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD20/CD8': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD4/CD45RA': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD4/CD8': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
#> Warning: Got the following message when running wilcox.test test for marker 'CD45RA/CD8': Sample1 vs Sample2
#>   cannot compute exact confidence intervals with ties
dca_markers
#> # A tibble: 10 × 15
#>      estimate data_type target  reference    n1    n2 statistic     p p_adj
#>         <dbl> <chr>     <chr>   <chr>     <int> <int>     <dbl> <dbl> <dbl>
#>  1  0         pearson_z Sample1 Sample2       4     4       8       1     1
#>  2  0         pearson_z Sample1 Sample2       5     5      12.5     1     1
#>  3 -0.0000165 pearson_z Sample1 Sample2       5     5      12.5     1     1
#>  4  0         pearson_z Sample1 Sample2       5     5      12.5     1     1
#>  5  0         pearson_z Sample1 Sample2       4     4       8       1     1
#>  6  0         pearson_z Sample1 Sample2       4     4       8       1     1
#>  7  0         pearson_z Sample1 Sample2       4     4       8       1     1
#>  8  0         pearson_z Sample1 Sample2       5     5      12.5     1     1
#>  9  0         pearson_z Sample1 Sample2       5     5      12.5     1     1
#> 10 -0.0000464 pearson_z Sample1 Sample2       5     5      12.5     1     1
#> # ℹ 6 more variables: conf.low <dbl>, conf.high <dbl>, method <chr>,
#> #   alternative <chr>, marker_1 <chr>, marker_2 <chr>
```
