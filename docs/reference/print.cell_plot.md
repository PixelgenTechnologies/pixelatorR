# Print a cell plot recipe

Builds and draws the recipe as a static ggplot, then returns that ggplot
invisibly. This matches `print.ggplot()`: typing a recipe at the prompt
renders the plot. Use
[`cell_plot_interactive()`](cell_plot_interactive.md) for an interactive
rgl or Plotly widget.

## Usage

``` r
# S3 method for class 'cell_plot'
print(x, ...)
```

## Arguments

- x:

  A `cell_plot` recipe.

- ...:

  Additional arguments. Currently not used.

## Value

A `ggplot` object, invisibly.
