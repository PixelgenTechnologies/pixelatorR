# Differential analysis of proximity scores

Runs differential analysis (Running Wilcoxon rank-sum test) on proximity
scores calculated from Proximity Network Assay (PNA) data generated with
the `Pixelator` data processing pipeline.

## Usage

``` r
DifferentialProximityAnalysis(object, ...)

# S3 method for class 'data.frame'
DifferentialProximityAnalysis(
  object,
  contrast_column,
  reference,
  targets = NULL,
  group_vars = NULL,
  proximity_metric = "join_count_z",
  metric_type = c("all", "self", "co"),
  backend = c("dplyr", "data.table"),
  min_n_obs = 0,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
)

# S3 method for class 'Seurat'
DifferentialProximityAnalysis(
  object,
  contrast_column,
  reference,
  targets = NULL,
  assay = NULL,
  group_vars = NULL,
  proximity_metric = "join_count_z",
  metric_type = c("all", "self", "co"),
  min_n_obs = 0,
  p_adjust_method = c("bonferroni", "holm", "hochberg", "hommel", "BH", "BY", "fdr"),
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  An object containing PNA proximity scores

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

- proximity_metric:

  The proximity metric to use. Any numeric data column in the proximity
  score table can be selected. The default is "pearson_z".

- metric_type:

  One of "all", "self" or "cross". If "all", all pairwise comparisons
  are considered. If "self", only protein pairs of the same type are
  considered. If "cross", protein pairs of different type are
  considered.

- backend:

  One of "dplyr" or "data.table". The latter requires the `dtplyr`
  package to be installed.

- min_n_obs:

  Minimum number of observations allowed in a group. Target groups with
  less observations than `min_n_obs` will be skipped.

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

If you are working with a `Seurat` object created with pixelatorR, the
proximity scores can be accessed with
[`ProximityScores`](ProximityScores.md).

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
Proximity Network Assay (PNA) data. `seurat_object` contains a
`meta.data` column called "sampleID" that holds information about what
samples the components originated from. This column could have three
sample IDs: "control", "stimulated1" and "stimulated2". In addition, we
have a column called "cell_type" that holds information about the cell
type identity of each component.

1.  If we want to compare the "stimulated1" group to the "control"
    group:

        dp_markers <- DifferentialProximityAnalysis(
           object = seurat_object,
           contrast_column = "sampleID",
           reference = "control",
           targets = "stimulated1"
        )

2.  If we want to compare the "stimulated1" and "stimulated2" groups to
    the "control" group:

        dp_markers <- DifferentialProximityAnalysis(
          object = seurat_object,
          contrast_column = "sampleID",
          reference = "control",
          targets = c("stimulated1", "stimulated2")
        )

3.  If we want to compare the "stimulated1" and "stimulated2" groups to
    the "control" group, and split the tests by cell type:

        dp_markers <- DifferentialProximityAnalysis(
           object = seurat_object,
           contrast_column = "sampleID",
           reference = "control",
           targets = c("stimulated1", "stimulated2"),
           group_vars = "cell_type"
        )

## See also

Other DA-methods: [`RunDAA()`](RunDAA.md), [`RunDCA()`](RunDCA.md),
[`RunDPA()`](RunDPA.md)

## Examples

``` r
# TODO: Update examples with real data
library(dplyr)
example_data <- tidyr::expand_grid(
  marker_1 = c("HLA-ABC", "B2M", "CD4", "CD8", "CD20", "CD19", "CD45", "CD43") %>%
    rep(each = 50),
  marker_2 = c("HLA-ABC", "B2M", "CD4", "CD8", "CD20", "CD19", "CD45", "CD43") %>%
    rep(each = 50)
) %>%
  mutate(
    join_count_z = rnorm(n(), sd = 10)
  )

example_data <- example_data %>%
  mutate(sampleID = "ctrl") %>%
  bind_rows(
    example_data %>% mutate(join_count_z = join_count_z + 1) %>%
      mutate(sampleID = "treatment")
  )

# Compute statistics
dp_results <- DifferentialProximityAnalysis(
  example_data,
  contrast_column = "sampleID",
  reference = "ctrl",
  proximity_metric = "join_count_z",
  metric_type = "self"
)
#> ℹ Computing Running Wilcoxon rank-sum test for each marker pair across the following comparisons:
#> 
#>   • treatment vs ctrl
```
