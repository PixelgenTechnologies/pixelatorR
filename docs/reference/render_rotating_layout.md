# Create a rotating 3D layout video

**\[experimental\]**

`render_rotating_layout` can be used to generate a rotating 3D scatter
plot from a tibble with node layout coordinates.
`render_rotating_layout` offers a number of customization options,
including the option to facet the plot by additional columns in the
tibble. See the sections below for more details on how to adjust the
appearance.

## Usage

``` r
render_rotating_layout(
  data,
  file,
  cell_col = NULL,
  marker_col = NULL,
  pt_opacity = 0.7,
  pt_size = 1,
  colors = RColorBrewer::brewer.pal(9, "Blues"),
  max_degree = 360,
  center_zero = FALSE,
  scale_layout = TRUE,
  frames = 500,
  pad = 0.1,
  show_first_frame = TRUE,
  width = 500,
  height = 500,
  res = 150,
  delay = 1/20,
  ggplot_theme = NULL,
  title = "",
  bg = "white",
  label_grid_axes = TRUE,
  margin_widths = c(0.1, 0.1),
  use_facet_grid = FALSE,
  flip = FALSE,
  graphics_use = c("base", "ggplot2"),
  boomerang = FALSE,
  cl = NULL,
  keep_frames = FALSE
)
```

## Arguments

- data:

  A tibble (`tbl_df`) with columns 'x', 'y', 'z', and 'node_val'. The
  'node_val' column can be either a numeric or a factor.

- file:

  A character string specifying the path to the output file. The video
  format will be determined by the file extension. If the extension is
  '.gif', the `gifski` R package will be used to render a gif. Other
  file formats (such as mp4, mkv, mov, or flv) will be rendered using
  the `av` R package. Make sure to use a file format supported by either
  `gifski` or `av`.

- cell_col, marker_col:

  A character string specifying columns to facet the plot by. `cell_col`
  should be a column with unique cell identifiers and `marker_col`
  should be a column with unique protein identifiers.

- pt_opacity:

  A numeric value between 0 and 1 specifying the opacity of the points.

- pt_size:

  A numeric value indicating the maximum size of the points.

- colors:

  A vector of valid colors. If 'node_val' is a factor, the length of
  'colors' should be equal to the number of unique levels. If 'node_val'
  is numeric, the colors will be used to create a gradient color scale.

- max_degree:

  A numeric value between 90 and 360. The maximum angle of rotation
  around the z axis. Default is 360, which corresponds to a full turn
  around the z axis.

- center_zero:

  A logical value indicating whether the color gradient should be
  centered around zero.

- scale_layout:

  A logical value indicating whether the layouts should be such that
  they all fit the same bounding box. This is typically necessary as
  layouts vary significantly in size.

- frames:

  A positive numeric value indicating the number of frames to render.
  More frames will result in a smoother animation but will increase the
  size of the output file and rendering time.

- pad:

  A numeric value between 0 and 1. The amount of padding to add to the
  axis range. Default is 0.1 which corresponds to an expansion of 10%
  along all three axes.

- show_first_frame:

  A logical value indicating whether the first frame should be displayed
  in the viewer before proceeding with the rendering. In interactive
  sessions, you will be prompted to enter 'Yes' to continue if you are
  satisfied with the first frame. Default is TRUE.

- width, height:

  An integer specifying the width and height of the video in pixels.
  Default is 500x500 which is suitable for a low resolution GIF file.

- res:

  A numeric value indicating the resolution of the png files in dpi. See
  [`png`](https://rdrr.io/r/grDevices/png.html) for more details.

- delay:

  A numeric value indicating the delay between frames in seconds.
  Default is 1/20, which corresponds to 20 frames per second. The value
  should be between 1/100 and 1.

- ggplot_theme:

  A ggplot2 theme object. Default is `NULL`, which corresponds to the
  default plot theme.

- title:

  A character string specifying the title of the plot.

- bg:

  A valid color string specifying the background color of the plot.

- label_grid_axes:

  A logical value indicating whether the rows and columns of the plot
  grid should be labeled.

- margin_widths:

  A numeric vector of length 2 specifying the width of the left and top
  margins relative to the width and height of the plot. Default is
  c(0.1, 0.1), meaning that the left and top margins will take up 10% of
  the width and height respectively.

- use_facet_grid:

  If set to TRUE, the plot will be faceted using
  [`facet_grid`](https://ggplot2.tidyverse.org/reference/facet_grid.html).
  Note that with this option, the node colors will be mapped to a single
  color scale. If these selected markers have very different dynamic
  ranges, it might be difficult to observe trends in lowly expressed
  markers.

- flip:

  A logical value indicating whether the plot should be flipped such
  that markers are arranged in columns and cells in rows.

- graphics_use:

  One of "ggplot2" or "base". The default, "base", is much faster but
  has fewer customization options. The "ggplot2" option comes with a
  number of customization options.

- boomerang:

  A logical value specifying whether the rotation should be reversed
  after reaching the maximum degree to create a "boomerang" effect. The
  frames are simply reversed and appended to the original frames. The
  resulting video will have twice as many frames as specified by the
  `frames` argument and keep the same frame rate. Note that if the
  output is a gif file, it will double in size.

- cl:

  A number of threads or a cluster object created by
  [`makeCluster`](https://rdrr.io/r/parallel/makeCluster.html). The
  default is `NULL`, which corresponds to a single-threaded operation
  (sequential processing).

- keep_frames:

  A logical value indicating whether the png files should be kept after
  rendering the video. This option is useful if you want to access to
  the individual frames later for further processing.

## Value

Exports an animation of a rotating 3D scatter plot.

## cells and markers

The `cell_col` and `marker_col` arguments can be used to facet the plot
by additional columns in the tibble. For instance, if the tibble
contains layout coordinates for multiple cells and markers, two
additional columns should specify these groupings. In the final plot,
each cell will plotted in its own coordinate system to make sure that
the layout fills the plot area. Similarly, each marker will share the
same color scale to make sure that the colors are comparable across
cells.

## graphic parameters

The appearance of the plot can be customized using a number of graphics
parameters. Increasing the number of `frames` while decreasing the delay
will create a smoother transition but will take much longer to render.
The `width` and `height` parameters should be increased manually to
accommodate more cells and markers or to increase resolution. The
`width` and `height` parameters should be adjusted to match the selected
dpi (`res`). Note that when these parameters are increased, the
rendering time can increase drastically. Unfortunately, drawing ggplot
objects to the graphics device is relatively slow so the rendering might
take several minutes. Using base R graphics instead of ggplot2 can speed
up the rendering significantly, but it will limit the customization
options.

Since the rendering is expensive, we might want to inspect the first
frame to make sure that the plot looks as expected. If the
`show_first_frame` argument is set to `TRUE`, which is the default, the
first frame will be displayed in the viewer before proceeding with the
rendering. In interactive sessions, you will be prompted to enter 'Yes'
to continue if you are satisfied with the first frame.

## parallel processing

To speed up the rendering, the `cl` argument can be used to specify a
number of threads. This feature only works on unix-based systems. This
can significantly reduce the rendering time; however, it might create
artefacts in the output video. The default is `NULL`, which corresponds
to a single-threaded operation (sequential processing).

## graphics_use

The `graphics_use` argument can be used to switch between the default
"base" and the "ggplot2" graphics system. The "base" option is much
faster but has fewer customization options. For instance, it is
currently not possible to get row and column labels with the "base"
graphics option.

## Examples

``` r
library(dplyr)
pxl_file <- minimal_pna_pxl_file()
se <- ReadPNA_Seurat(pxl_file)
se <- se %>%
  LoadCellGraphs(add_layouts = TRUE)

# Create a gif from a 3D layout
cg <- CellGraphs(se)[[3]]
df <- cg@layout$wpmds_3d %>%
  mutate(node_val = cg@counts[, "CD3e"])
temp_gif <- fs::file_temp(ext = ".gif")
render_rotating_layout(
  data = df,
  file = temp_gif,
  colors = c("lightgrey", "red"),
  pt_size = 0.5,
  max_degree = 90,
  frames = 20,
  delay = 1 / 5,
  bg = "transparent",
  label_grid_axes = FALSE,
  show_first_frame = FALSE
)
#> ! The following parameters are not supported with graphics_use = 'base':
#> ggplot_theme, label_grid_axes, title
#> Titles, text annotations and color bar will be missing from the plot.
#> 
#> ── Rendering frames... 
#> 
#> ── Creating video... 
#> ✔ The video has been saved to:
#> C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf45b033fd3.gif
magick::image_read(temp_gif)
#> # A tibble: 20 × 7
#>    format width height colorspace matte filesize density
#>    <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
#>  1 GIF      500    500 sRGB       TRUE         0 72x72  
#>  2 GIF      500    500 sRGB       TRUE         0 72x72  
#>  3 GIF      500    500 sRGB       TRUE         0 72x72  
#>  4 GIF      500    500 sRGB       TRUE         0 72x72  
#>  5 GIF      500    500 sRGB       TRUE         0 72x72  
#>  6 GIF      500    500 sRGB       TRUE         0 72x72  
#>  7 GIF      500    500 sRGB       TRUE         0 72x72  
#>  8 GIF      500    500 sRGB       TRUE         0 72x72  
#>  9 GIF      500    500 sRGB       TRUE         0 72x72  
#> 10 GIF      500    500 sRGB       TRUE         0 72x72  
#> 11 GIF      500    500 sRGB       TRUE         0 72x72  
#> 12 GIF      500    500 sRGB       TRUE         0 72x72  
#> 13 GIF      500    500 sRGB       TRUE         0 72x72  
#> 14 GIF      500    500 sRGB       TRUE         0 72x72  
#> 15 GIF      500    500 sRGB       TRUE         0 72x72  
#> 16 GIF      500    500 sRGB       TRUE         0 72x72  
#> 17 GIF      500    500 sRGB       TRUE         0 72x72  
#> 18 GIF      500    500 sRGB       TRUE         0 72x72  
#> 19 GIF      500    500 sRGB       TRUE         0 72x72  
#> 20 GIF      500    500 sRGB       TRUE         0 72x72  

if (FALSE) { # \dontrun{
# Include multiple facets
markers <- paste0("Marker", 1:3)
df <- lapply(1:2, function(i) {
  set.seed(i)
  cg <- simulate_bipartite_graph(n_nodes = 5e3, epsilon = 8) %>%
    add_binary_marker_counts_pol(marker_polarize = "Marker1", epsilon = 12)
  xyz <- layout_with_weighted_pmds(cg@cellgraph, dim = 3) %>%
    as_tibble(.name_repair = ~ c("x", "y", "z"))
  lapply(colnames(cg@counts)[1:2], function(m) {
    xyz %>%
      mutate(node_val = cg@counts[, m] %>% as.character()) %>%
      mutate(marker = m)
  }) %>%
    bind_rows() %>%
    mutate(cell = paste0("Cell", i))
}) %>% bind_rows()

temp_gif <- fs::file_temp(ext = ".gif")
render_rotating_layout(
  data = df,
  pt_size = 0.5,
  width = 740,
  height = 650,
  cell_col = "cell",
  marker_col = "marker",
  colors = viridis::viridis(n = 2),
  file = temp_gif,
  max_degree = 90,
  frames = 100,
  delay = 1 / 20,
  res = 150,
  bg = "grey",
  show_first_frame = FALSE,
  cl = 9
)
magick::image_read(temp_gif)

# Use base R graphics
temp_gif <- fs::file_temp(ext = ".gif")
render_rotating_layout(
  data = df,
  pt_size = 0.5,
  pt_opacity = 1,
  width = 740,
  height = 650,
  cell_col = "cell",
  marker_col = "marker",
  colors = viridis::viridis(n = 2),
  file = temp_gif,
  max_degree = 360,
  frames = 500,
  delay = 1 / 20,
  res = 150,
  bg = "grey",
  show_first_frame = FALSE,
  graphics_use = "base",
  cl = 9
)
magick::image_read(temp_gif)
} # }
```
