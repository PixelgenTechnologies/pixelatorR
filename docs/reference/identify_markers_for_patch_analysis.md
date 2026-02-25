# Identify population-specific markers for patch detection

Patch detection can be applied to identify patches of one cell type on
another. The prerequisite for this analysis is that the patch-specific
markers are known. Patch detection can be improved further by also
leveraging receiver-specific markers. This function identifies
population-specific markers for the patch and receiver populations,
which can then be used for patch detection.

## Usage

``` r
identify_markers_for_patch_analysis(
  object,
  group_by,
  receiver_population,
  target_population,
  abundance_difference = 10,
  min_freq = 0.01,
  show_plot = TRUE,
  seed = 123
)
```

## Arguments

- object:

  A `Seurat` object

- group_by:

  A string specifying the metadata column to group by

- receiver_population:

  A string specifying the receiver population name present in the
  `group_by` column

- target_population:

  A string specifying the target population name present in the
  `group_by` column

- abundance_difference:

  A numeric value specifying how many times higher the abundance of a
  marker should be in one population relative to the other. This is only
  used to label the markers in the output.

- min_freq:

  A numeric value specifying the minimum frequency of a protein to be
  labeled in the output table.

- show_plot:

  Logical, whether to show a plot summarizing the results.

- seed:

  An integer seed for reproducibility

## Value

A `tbl_df` with the following columns:

- `marker`: the name of the protein.

- `receiver_unmixed_freq`: the estimated proportion in the receiver
  population after unmixing.

- `target_unmixed_freq`: the estimated proportion in the target
  population after unmixing.

- `receiver_freq`: the proportion in the receiver population.

- `target_freq`: the proportion in the target population.

- `label`: a label indicating whether the protein is a marker for the
  receiver or target population. NA values indicate that the protein is
  unspecific.

## Details

Patch detection is sensitive to the selection of markers and therefore
requires careful selection. The best markers are both high-abundant and
specific to a cell population.

The method requires a `Seurat` object with a metadata column containing
the population information, e.g. a column with cell type labels. We then
need to specify the receiver and target populations, where the receiver
population represent the cells on which the patches are expected to be
found and the target population represents the cell type from which the
patches originated.

As a practical example, let's say that our data represent co-cultured T
and B cells, and we anticipate that the T cells have patches of B cells
on them. Now we face a challenge because the T cell population contains
a lot of B cell markers, making it harder to determine what markers are
T-cell specific. In other words, the abundance data is mixed.

This method attempts to unmix the abundance data using matrix
factorization, estimating the composition of the pure receiver and
target populations. The unmixed abundance profiles are then used to
label population-specific markers based on the difference in abundance
and minimum frequency. The results are summarized in a table, and an
optional plot is drawn to help interpret the results.
