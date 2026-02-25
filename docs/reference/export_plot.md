# Save a ggplot2 plot to one or more file formats with options for size and directory creation.

This function is a wrapper around
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html) that
allows for saving the same plot in multiple formats and handles
directory creation if needed. It also includes an option to skip saving
if the file already exists, based on a global option.

## Usage

``` r
export_plot(
  filename,
  plot = ggplot2::last_plot(),
  width = 10,
  height = 10,
  create_dir = TRUE,
  file_formats = getOption("export_plot.file_formats", c("png", "pdf")),
  overwrite = getOption("export_plot.overwrite", FALSE)
)
```

## Arguments

- filename:

  Character. The path and base file name for saving the plot. If a file
  extension is included, it will override `file_formats` and save only
  in that format.

- plot:

  ggplot object. The plot to save. Defaults to the last plot produced
  ([`last_plot()`](https://ggplot2.tidyverse.org/reference/get_last_plot.html)).

- width:

  Numeric. The width of the plot in inches. Default is 10.

- height:

  Numeric. The height of the plot in inches. Default is 10.

- create_dir:

  Logical. Whether to create the directory if it does not exist. Default
  is TRUE.

- file_formats:

  Character vector. File formats to save to (e.g., "png", "pdf").
  Default is c("png", "pdf").

- overwrite:

  Logical. Whether to overwrite existing files. Default is taken from
  the global option "export_plot.overwrite" (default FALSE).

## Value

Invisibly returns NULL. Used for its side effect of writing plot files.

## Examples

``` r
if (FALSE) { # \dontrun{
export_plot("results/myplot", plot = myplot, width = 8, height = 6)
} # }
# Change overwrite behavior globally
options(export_plot.overwrite = TRUE)

# Now existing files will be overwritten without needing to specify overwrite in the function call
if (FALSE) { # \dontrun{
export_plot("results/myplot", plot = myplot, width = 8, height = 6)
} # }
```
