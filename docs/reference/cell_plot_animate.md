# Render a rotating cell plot animation

Builds a [`cell_plot()`](cell_plot.md) recipe once, draws one frame per
rotation angle, and encodes a GIF or video. Rotation geometry comes from
[`cell_coord_rotate()`](cell_coord_rotate.md). File type, size,
resolution, frame rate, and the frame backend belong here.

## Usage

``` r
cell_plot_animate(
  object,
  file,
  frames = 500,
  width = 500,
  height = 500,
  res = 150,
  fps = 20,
  frame_backend = c("base", "ggplot2"),
  workers = 1L
)
```

## Arguments

- object:

  A `cell_plot` recipe that includes
  [`cell_coord_rotate()`](cell_coord_rotate.md).

- file:

  Output path. The extension selects the encoder.

- frames:

  Positive whole number of encoded frames.

- width, height:

  Output size in pixels.

- res:

  PNG resolution in pixels per inch.

- fps:

  Encoded frames per second.

- frame_backend:

  Device used to draw each frame. `"base"` is faster. `"ggplot2"`
  matches the static renderer more closely.

- workers:

  Positive whole number of parallel workers. `1` renders sequentially.

## Value

The output path, invisibly.

## Details

GIF output uses gifski. Any other extension is encoded with av. The
parent directory of `file` must already exist. An existing file is
overwritten.

## See also

[`cell_plot()`](cell_plot.md),
[`cell_coord_rotate()`](cell_coord_rotate.md),
[`cell_illuminate()`](cell_illuminate.md)

## Examples

``` r
se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#> ✔ Created a <Seurat> object with 5 cells and 158 targeted surface proteins
se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
  ComputeLayout(layout_method = "spectral")
#> ℹ Computing layouts for 1 graphs

cell_graph <- CellGraphs(se)[[4]]

layout_data <- FetchLayoutData(cell_graph, layout_method = "spectral_3d") |>
  # Downsample to speed up tests
  dplyr::slice_sample(n = 5000)
file <- tempfile(fileext = ".gif")
cell_plot(layout_data) |>
  cell_coord_rotate(axis = "y") |>
  cell_plot_animate(file, frames = 2, width = 160, height = 160, res = 72)
#> ℹ Encoding /tmp/RtmpjKKHFf/file4d6b27fa317b.gif
#> ✔ Encoding /tmp/RtmpjKKHFf/file4d6b27fa317b.gif [18ms]
#> 
```
