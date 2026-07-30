library(dplyr)
prox <- ReadPNA_proximity(minimal_pna_pxl_file())
prox_summarized <- prox %>%
  slice_sample(n = 1e4) %>%
  group_by(marker_1, marker_2) %>%
  summarize(mean_log2_ratio = mean(log2_ratio), .groups = "drop") %>%
  mutate(test1 = marker_1, test2 = marker_2, estimate = mean_log2_ratio) %>%
  ungroup()

test_that("ColocalizationHeatmap works as expected", {
  # Default method
  expect_no_error(ColocalizationHeatmap(prox_summarized))

  # Dots method
  expect_no_error(p <- ColocalizationHeatmap(prox_summarized, type = "dots", size_col = "mean_log2_ratio"))
  expect_true(is.factor(p$data$marker_1))
  expect_true(is.factor(p$data$marker_2))
  expect_equal(p$data$marker_1 %>% levels(), p$data$marker_2 %>% levels())

  # Custom marker_1 and marker_2 columns
  expect_no_error(ColocalizationHeatmap(prox_summarized, marker1_col = "test1", marker2_col = "test2"))

  # Custom value column
  expect_no_error(ColocalizationHeatmap(prox_summarized, value_col = "mean_log2_ratio"))

  # Symmetrise FALSE
  expect_no_error(p <- ColocalizationHeatmap(prox_summarized, symmetrise = FALSE, type = "dots", size_col = "mean_log2_ratio"))
  expect_true(is.factor(p$data$marker_1))
  expect_true(is.factor(p$data$marker_2))
  expect_true(!identical(p$data$marker_1 %>% levels(), p$data$marker_2 %>% levels()))

  # highlight_pairs = NULL is a no-op (including default)
  expect_no_error(ColocalizationHeatmap(prox_summarized))
  expect_no_error(ColocalizationHeatmap(prox_summarized, highlight_pairs = NULL))

  # Outline selected pairs (dots + tiles)
  markers <- unique(c(prox_summarized$marker_1, prox_summarized$marker_2))
  highlight_pairs <- tibble(
    marker_1 = markers[1],
    marker_2 = markers[min(2, length(markers))]
  )
  expect_no_error(
    p <- ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_pairs,
      highlight_colors = "black"
    )
  )
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomTile"), logical(1))))
  tile_layer <- which(vapply(
    p$layers,
    function(l) inherits(l$geom, "GeomTile"),
    logical(1)
  ))[1]
  tile_data <- ggplot_build(p)$data[[tile_layer]]
  expect_equal(unique(round(tile_data$xmax - tile_data$xmin, 10)), 0.9)
  expect_equal(unique(round(tile_data$ymax - tile_data$ymin, 10)), 0.9)

  expect_no_error(
    p_shrunk <- ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_pairs,
      highlight_shrink = 0.25
    )
  )
  tile_layer <- which(vapply(
    p_shrunk$layers,
    function(l) inherits(l$geom, "GeomTile"),
    logical(1)
  ))[1]
  tile_data <- ggplot_build(p_shrunk)$data[[tile_layer]]
  expect_equal(unique(round(tile_data$xmax - tile_data$xmin, 10)), 0.75)
  expect_equal(unique(round(tile_data$ymax - tile_data$ymin, 10)), 0.75)

  # With symmetrise = FALSE, dots outline only the given orientation
  expect_no_error(
    p_asym <- ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      symmetrise = FALSE,
      highlight_pairs = highlight_pairs
    )
  )
  tile_layer <- which(vapply(p_asym$layers, function(l) inherits(l$geom, "GeomTile"), logical(1)))[1]
  expect_equal(nrow(p_asym$layers[[tile_layer]]$data), 1)
  expect_no_error(
    ColocalizationHeatmap(
      prox_summarized,
      type = "tiles",
      highlight_pairs = highlight_pairs,
      highlight_colors = "#C0392B"
    )
  )
  expect_no_error(
    ColocalizationHeatmap(
      prox_summarized,
      type = "tiles",
      symmetrise = FALSE,
      highlight_pairs = highlight_pairs
    )
  )

  # Discrete border colours from highlight_color_col (dots)
  highlight_disc <- highlight_pairs %>%
    mutate(database = "string")
  expect_no_error(
    p_disc <- ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_disc,
      highlight_color_col = "database",
      highlight_colors = c(string = "#C0392B", alphafold = "#6C3483")
    )
  )
  expect_true(any(vapply(p_disc$scales$scales, function(s) {
    inherits(s, "ScaleDiscrete") && "colour" %in% s$aesthetics
  }, logical(1))))

  # Continuous border colours from highlight_color_col (dots)
  highlight_cont <- highlight_pairs %>%
    mutate(score = 0.8)
  expect_no_error(
    p_cont <- ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_cont,
      highlight_color_col = "score",
      highlight_colors = c("#FEE5D9", "#A50F15")
    )
  )
  expect_true(any(vapply(p_cont$scales$scales, function(s) {
    inherits(s, "ScaleContinuous") && "colour" %in% s$aesthetics
  }, logical(1))))

  # Factor marker columns must still match character highlight_pairs
  prox_factor <- prox_summarized %>%
    mutate(
      marker_1 = factor(marker_1),
      marker_2 = factor(marker_2)
    )
  expect_no_error(
    ColocalizationHeatmap(
      prox_factor,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_pairs
    )
  )
})

test_that("ColocalizationHeatmap fails with invalid input", {
  expect_error(ColocalizationHeatmap("Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, marker1_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, marker2_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, value_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, type = "dots", size_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, size_col_transform = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, size_range = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, colors = FALSE))
  expect_error(ColocalizationHeatmap(prox_summarized, cluster_rows = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, cluster_cols = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, type = "Invalid"))
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      highlight_pairs = tibble(a = 1, b = 2)
    )
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      highlight_shrink = -0.1
    ),
    "at least 0 and less than 1"
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      highlight_shrink = 1
    ),
    "at least 0 and less than 1"
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      highlight_shrink = NA_real_
    ),
    "at least 0 and less than 1"
  )

  markers <- unique(c(prox_summarized$marker_1, prox_summarized$marker_2))
  highlight_pairs <- tibble(
    marker_1 = markers[1],
    marker_2 = markers[min(2, length(markers))],
    database = "string",
    score = 0.5
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      highlight_colors = c("black", "red")
    ),
    "highlight_color_col"
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      type = "tiles",
      highlight_pairs = highlight_pairs,
      highlight_color_col = "database",
      highlight_colors = c(string = "black")
    ),
    "dots"
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_pairs,
      highlight_color_col = "database",
      highlight_colors = c(alphafold = "#6C3483")
    ),
    "missing named colours"
  )
  expect_error(
    ColocalizationHeatmap(
      prox_summarized,
      type = "dots",
      size_col = "mean_log2_ratio",
      highlight_pairs = highlight_pairs,
      highlight_color_col = "score",
      highlight_colors = "red"
    ),
    "scale_color_gradientn"
  )
})
