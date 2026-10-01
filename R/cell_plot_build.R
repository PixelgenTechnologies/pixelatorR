#' Build a cell plot recipe
#'
#' Resolves mappings, constants, scales, illumination, and display defaults
#' without drawing a plot. The returned object has the same top-level fields as
#' the input recipe, with data ordered by the `arrange` mapping. Rows with
#' missing values (`NA`) in mapped `size` or `alpha` columns are removed with a
#' warning before color illumination is calculated. Color, size, and alpha
#' resolve to one value each when the aesthetic is constant or unmapped. Named
#' categorical size and alpha scales assign one resolved value per level.
#' Illumination is calculated after color-scale training and stored separately
#' from the unilluminated resolved colors used by legends. Other
#' projection-specific transformations, such as apparent size from depth, are
#' applied by renderers rather than during the build.
#'
#' @param object A `cell_plot` recipe.
#'
#' @return A `cell_plot_built` object containing resolved plot data and
#' rendering instructions.
#'
#' @examples
#' se <- ReadPNA_Seurat(minimal_pna_pxl_file())
#' se <- LoadCellGraphs(se, cells = colnames(se)[4], verbose = FALSE) |>
#'   ComputeLayout(layout_method = "spectral")
#'
#' cell_graph <- CellGraphs(se)[[4]]
#'
#' layout_data <- FetchLayoutData(cell_graph, vars = "CD82", layout_method = "spectral_3d") |>
#'   # Downsample to speed up tests
#'   dplyr::slice_sample(n = 5000)
#'
#' built <- cell_plot(layout_data, color = "CD82") |>
#'   build_cell_plot()
#'
#' @export
build_cell_plot <- function(object) {
  .validate_cell_plot(object)
  expect_scales()

  mapping <- object$mapping
  data <- object$data

  if (!is.null(mapping$arrange)) {
    data <- data[
      order(data[[mapping$arrange]], na.last = TRUE),
    ]
  }

  size_values <- if (is.null(mapping$size)) {
    NULL
  } else {
    data[[mapping$size]]
  }
  size_scale <- .build_cell_numeric_scale(
    values = size_values,
    specification = object$size,
    output_name = "sizes",
    default_output = c(2, 6),
    constant = object$constant$size %||% 1
  )

  alpha_values <- if (is.null(mapping$alpha)) {
    NULL
  } else {
    data[[mapping$alpha]]
  }
  alpha_scale <- .build_cell_numeric_scale(
    values = alpha_values,
    specification = object$alpha,
    output_name = "alphas",
    default_output = c(0.2, 1),
    constant = object$constant$alpha %||% 1
  )

  dropped <- .cell_plot_drop_missing_size_alpha(
    data = data,
    mapping = mapping,
    size_scale = size_scale,
    alpha_scale = alpha_scale
  )
  data <- dropped$data
  size_scale <- dropped$size_scale
  alpha_scale <- dropped$alpha_scale

  color_values <- if (is.null(mapping$color)) {
    NULL
  } else {
    data[[mapping$color]]
  }
  if (!is.null(color_values)) {
    has_usable_color <- if (is.numeric(color_values)) {
      any(is.finite(color_values))
    } else {
      any(!is.na(color_values))
    }
    if (!has_usable_color) {
      cli::cli_abort(
        c(
          "x" = paste(
            "Mapped {.arg color} column {.str {mapping$color}} has no",
            "usable values after rows with missing {.field size} or",
            "{.field alpha} were removed."
          )
        )
      )
    }
  }
  color_scale <- .build_cell_color_scale(
    values = color_values,
    specification = object$color,
    constant = object$constant$color
  )
  color_scale["illuminated"] <- list(NULL)
  if (!is.null(object$illuminate)) {
    layout <- tibble::tibble(
      x = data[[mapping$x]],
      y = data[[mapping$y]],
      z = data[[mapping$z]]
    )
    color_scale$illuminated <- .cell_illuminated_colors(
      layout = layout,
      colors = color_scale$resolved,
      mapped_values = color_scale$values,
      panel_rows = .cell_plot_panel_rows(data, object$grid),
      specification = object$illuminate,
      illumination_mask = if (is.null(mapping$illumination_mask)) {
        NULL
      } else {
        data[[mapping$illumination_mask]]
      }
    )
  }

  built <- structure(
    list(
      data = data,
      mapping = mapping,
      constant = object$constant,
      grid = object$grid,
      color = color_scale,
      size = size_scale,
      depth = object$depth,
      alpha = alpha_scale,
      theme = object$theme %||% .cell_plot_default_theme,
      annotation = object$annotation %||% .cell_plot_default_annotation,
      illuminate = object$illuminate,
      coord = object$coord
    ),
    class = "cell_plot_built"
  )

  return(built)
}

#' Drop rows with missing mapped size or alpha
#'
#' Removes rows whose resolved size or alpha is missing, keeping the mapped
#' values and resolved output in sync with `data`. Call this before color
#' illumination so dropped rows do not affect lighting.
#'
#' @param data Arranged cell plot data.
#' @param mapping Column mappings.
#' @param size_scale,alpha_scale Resolved numeric scales.
#'
#' @return A list with filtered `data`, `size_scale`, and `alpha_scale`.
#'
#' @noRd
.cell_plot_drop_missing_size_alpha <- function(
  data,
  mapping,
  size_scale,
  alpha_scale
) {
  size_drop <- if (!is.null(mapping$size)) is.na(size_scale$resolved) else FALSE
  alpha_drop <- if (!is.null(mapping$alpha)) is.na(alpha_scale$resolved) else FALSE
  drop_rows <- size_drop | alpha_drop

  if (!any(drop_rows)) {
    return(list(
      data = data,
      size_scale = size_scale,
      alpha_scale = alpha_scale
    ))
  }
  if (all(drop_rows)) {
    cli::cli_abort(
      "All rows were dropped due to missing values in {.field size} and/or {.field alpha}."
    )
  }

  cli::cli_warn(
    "Removed {sum(drop_rows)} row{?s} containing missing values in {.field size} and/or {.field alpha}."
  )

  keep_rows <- !drop_rows
  data <- data[keep_rows, , drop = FALSE]
  if (!is.null(mapping$size)) {
    size_scale <- .subset_cell_numeric_scale(size_scale, keep_rows)
  }
  if (!is.null(mapping$alpha)) {
    alpha_scale <- .subset_cell_numeric_scale(alpha_scale, keep_rows)
  }

  return(list(
    data = data,
    size_scale = size_scale,
    alpha_scale = alpha_scale
  ))
}

#' Subset resolved size or alpha scale vectors
#'
#' @param scale A numeric scale list.
#' @param keep_rows Logical rows to retain.
#'
#' @return The scale with `values` and `resolved` subset.
#'
#' @noRd
.subset_cell_numeric_scale <- function(scale, keep_rows) {
  scale$resolved <- scale$resolved[keep_rows]
  if (!is.null(scale$values)) {
    scale$values <- scale$values[keep_rows]
  }
  return(scale)
}

#' Find the data rows belonging to each cell plot panel
#'
#' Groups rows using the panel columns recorded by [cell_grid()]. Without a
#' grid, all rows belong to one panel.
#'
#' @param data Cell plot data.
#' @param grid Optional panel grid specification.
#'
#' @return A list of integer row indices, one element per non-empty panel.
#'
#' @noRd
.cell_plot_panel_rows <- function(data, grid) {
  panel_columns <- unname(unlist(grid[c("rows", "cols")], use.names = FALSE))
  if (length(panel_columns) == 0) {
    return(list(seq_len(nrow(data))))
  }

  rows <- data |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(panel_columns)),
      .drop = TRUE
    ) |>
    dplyr::group_rows()
  return(rows)
}

#' Apply geometry-based illumination to point colors
#'
#' Computes illumination independently within each non-empty plot panel. The
#' caller supplies coordinates so the same operation can be applied to original
#' or transformed layouts without retraining the color scale.
#'
#' @param layout A tibble with numeric `x`, `y`, and `z` columns.
#' @param colors Unilluminated hexadecimal colors, one per row or length one.
#' @param mapped_values Original mapped color values, or `NULL`.
#' @param panel_rows A list of integer row indices for non-empty panels.
#' @param specification An illumination specification.
#' @param illumination_mask Optional logical vector selecting rows whose colors
#' are illuminated.
#'
#' @return An illuminated hexadecimal color for every row in `layout`.
#'
#' @noRd
.cell_illuminated_colors <- function(
  layout,
  colors,
  mapped_values,
  panel_rows,
  specification,
  illumination_mask = NULL
) {
  illuminated_colors <- rep_len(colors, nrow(layout))

  for (rows in panel_rows) {
    if (length(rows) < 2) {
      next
    }

    illumination <- .cell_heuristic_illumination(
      layout = layout[rows, , drop = FALSE],
      clamp_quantiles = specification$clamp_quantiles,
      directional_light_weight = specification$directional_light_weight,
      volume_shading_weight = specification$volume_shading_weight,
      ambient_occlusion_weight = specification$ambient_occlusion_weight,
      ambient_occlusion_k = min(
        specification$ambient_occlusion_k,
        length(rows) - 1L
      ),
      normalize_weights = TRUE,
      light_direction = specification$light_direction
    )

    illumination <- .cell_normalize_illumination(illumination)

    apply_illumination <- rep(TRUE, length(rows))
    if (!is.null(mapped_values)) {
      apply_illumination <- !is.na(mapped_values[rows])
    }
    if (!is.null(illumination_mask)) {
      apply_illumination <- apply_illumination & illumination_mask[rows]
    }
    color_rows <- rows[apply_illumination]
    color_mask <- illumination[apply_illumination]
    if (length(color_rows) == 0) {
      next
    }

    if (is.null(specification$shadow_colors)) {
      illuminated_colors[color_rows] <- .cell_apply_hsv_illumination(
        illuminated_colors[color_rows],
        illumination = color_mask,
        ambient_intensity = specification$ambient_intensity,
        saturation_boost = specification$saturation_boost
      )
    } else {
      shadow_color <- scales::col_numeric(
        palette = specification$shadow_colors,
        domain = c(0, 1),
        na.color = "transparent"
      )(color_mask)
      illuminated_colors[color_rows] <- .cell_apply_palette_illumination(
        illuminated_colors[color_rows],
        shadow_color,
        illumination = color_mask,
        ambient_intensity = specification$ambient_intensity
      )
    }
  }

  return(illuminated_colors)
}

#' Select rendered colors from a built cell plot
#'
#' Returns baked illumination when available and otherwise returns the
#' unilluminated colors resolved by the color scale.
#'
#' @param object A `cell_plot_built` object.
#'
#' @return A hexadecimal color vector.
#'
#' @noRd
.cell_plot_rendered_colors <- function(object) {
  return(object$color$illuminated %||% object$color$resolved)
}

#' Apparent diameters from projected depth
#'
#' Scales point **diameter** with inverse camera distance so a point twice as
#' far from the camera is drawn at half the diameter. Larger `z` is closer.
#' Camera distance is `focal_distance + (z_max - z)`. Diameters are then
#' multiplied so the point at `mean(z)` keeps `base_size`. True 3D renderers
#' should not apply this factor. ggplot2 applies it only when size is unmapped.
#'
#' @param z Numeric depth coordinates.
#' @param base_size Constant relative size of the mean-depth point.
#' @param focal_distance Positive focal distance in the units of `z`.
#'
#' @return A numeric vector of apparent diameters, the same length as `z`.
#'
#' @noRd
.cell_depth_apparent_size <- function(
  z,
  base_size,
  focal_distance = .cell_plot_default_focal_distance
) {
  base_size <- base_size[[1]]
  if (length(z) == 0L) {
    return(z)
  }

  z_range <- range(z, na.rm = TRUE)
  if (!is.finite(diff(z_range)) || diff(z_range) == 0) {
    return(rep(base_size, length(z)))
  }

  z_max <- z_range[2]
  s <- focal_distance / (focal_distance + (z_max - z))
  s_mean <- focal_distance /
    (focal_distance + (z_max - mean(z, na.rm = TRUE)))
  return(base_size * s / s_mean)
}

#' Resolve relative point sizes for a projected cell plot
#'
#' Applies apparent depth sizing when a depth mapping is active and node size
#' is not mapped. Animation frames remap depth to their rotated z coordinate
#' before calling this helper.
#'
#' @param object A `cell_plot_built` object.
#'
#' @return Relative point diameters for the projected plot.
#'
#' @noRd
.cell_plot_projected_sizes <- function(object) {
  relative_size <- object$size$resolved
  if (is.null(object$mapping$depth) || !is.null(object$mapping$size)) {
    return(relative_size)
  }

  relative_size <- .cell_depth_apparent_size(
    object$data[[object$mapping$depth]],
    base_size = relative_size,
    focal_distance = object$depth$focal_distance %||%
      .cell_plot_default_focal_distance
  )
  return(relative_size)
}

#' Convert relative node sizes to ggplot2 millimetres
#'
#' One relative unit is one millimetre in ggplot2 point sizes.
#'
#' @param sizes Relative node sizes.
#'
#' @return ggplot2 point sizes in millimetres.
#'
#' @noRd
.cell_relative_size_to_ggplot <- function(sizes) {
  return(sizes)
}

#' Unique facet levels including unused factor levels
#'
#' @param values A facet column.
#'
#' @return A list of level values in display order.
#'
#' @noRd
.cell_plot_facet_levels <- function(values) {
  if (is.factor(values)) {
    levels <- as.list(levels(values))
    if (anyNA(values)) {
      levels <- c(levels, list(NA))
    }
    return(levels)
  }
  unique_values <- sort(unique(values), na.last = TRUE)
  return(as.list(unique_values))
}

#' Match rows belonging to one facet level
#'
#' @param values A facet column.
#' @param level The level to match, which may be `NA`.
#'
#' @return A logical vector.
#'
#' @noRd
.cell_plot_facet_match <- function(values, level) {
  if (is.na(level)) {
    return(is.na(values))
  }
  return(!is.na(values) & as.character(values) == as.character(level))
}

#' Convert a facet level to strip text
#'
#' @param level A facet level, which may be `NA`.
#'
#' @return A character label.
#'
#' @noRd
.cell_plot_facet_label <- function(level) {
  if (is.na(level)) {
    return("NA")
  }
  return(as.character(level))
}

#' Add proportional padding to a cell-plot coordinate range
#'
#' @param values Numeric coordinate values.
#'
#' @return A length-two numeric range with five percent padding.
#'
#' @noRd
.cell_plot_padded_range <- function(values) {
  limits <- range(values, na.rm = TRUE)
  pad <- if (diff(limits) > 0) {
    diff(limits) * 0.05
  } else {
    max(abs(limits[1]), 1) * 0.05
  }
  return(limits + c(-pad, pad))
}

#' Convert relative node sizes to screen pixel diameters
#'
#' Assumes 96 pixels per inch, matching Plotly CSS pixels and rgl point sizes.
#'
#' @param sizes Relative node sizes.
#'
#' @return Marker diameters in pixels.
#'
#' @noRd
.cell_relative_size_to_pixels <- function(sizes) {
  return(sizes * (96 / 25.4))
}

#' Convert relative node sizes to base R `cex`
#'
#' One relative unit is one millimetre. `cex = 1` is the default point size
#' (`ps`), which is 12 points unless a renderer supplies another value.
#'
#' @param sizes Relative node sizes.
#' @param ps Base point size corresponding to `cex = 1`.
#'
#' @return A `cex` multiplier.
#'
#' @noRd
.cell_relative_size_to_cex <- function(sizes, ps = 12) {
  millimetres_per_cex <- ps * 25.4 / 72
  return(sizes / millimetres_per_cex)
}

#' Resolve a requested color scale type
#'
#' @param values Mapped color values.
#' @param type `"auto"`, `"continuous"`, or `"categorical"`.
#'
#' @return `"continuous"` or `"categorical"`.
#'
#' @noRd
.cell_color_scale_type <- function(values, type = "auto") {
  if (identical(type, "continuous") || identical(type, "categorical")) {
    if (identical(type, "continuous") && !is.numeric(values)) {
      cli::cli_abort(
        c("x" = "Continuous color scales require a numeric or integer column.")
      )
    }
    return(type)
  }

  if (is.numeric(values)) {
    return("continuous")
  }
  return("categorical")
}

#' Build a cell color scale
#'
#' Resolves portable hex colors in `resolved` for backends that cannot train a
#' legend scale. ggplot retrains color from the mapped column and `colors` /
#' `limits` so it can draw a legend. Training always uses the full data,
#' including every panel. Without a mapped column, every point takes the
#' constant color and the scale carries no legend metadata.
#'
#' @param values Mapped color values or `NULL`.
#' @param specification A color specification or `NULL`.
#' @param constant A constant color used when no column is mapped.
#'
#' @return Resolved colors and scale metadata.
#'
#' @noRd
.build_cell_color_scale <- function(values, specification, constant = NULL) {
  default_colors <- c("lightgrey", "mistyrose", "red", "darkred")
  if (is.null(values)) {
    resolved <- .cell_colors_to_hex(constant %||% "gray90")
    return(list(
      colors = resolved,
      limits = NULL,
      na_color = NULL,
      type = NULL,
      values = NULL,
      resolved = resolved
    ))
  }

  specification <- specification %||% list()
  colors <- specification$colors %||% default_colors
  limits <- specification$limits
  na_color <- specification$na_color %||% "grey50"
  scale_type <- .cell_color_scale_type(values, specification$type %||% "auto")

  trained <- .resolve_cell_color_values(
    values = values,
    colors = colors,
    limits = limits,
    na_color = na_color,
    scale_type = scale_type
  )
  return(list(
    colors = trained$colors,
    limits = trained$limits,
    na_color = trained$na_color,
    type = trained$type,
    values = values,
    resolved = trained$resolved
  ))
}

#' Resolve colors for one set of mapped values
#'
#' @param values Mapped color values for one scale domain.
#' @param colors Palette colors.
#' @param limits Optional explicit limits.
#' @param na_color Missing-value color.
#' @param scale_type `"continuous"` or `"categorical"`.
#'
#' @return Palette, limits, type, and resolved colors.
#'
#' @noRd
.resolve_cell_color_values <- function(
  values,
  colors,
  limits,
  na_color,
  scale_type
) {
  if (scale_type == "continuous") {
    if (!any(is.finite(values))) {
      cli::cli_abort(
        c(
          "x" = "A continuous color scale needs at least one finite mapped value."
        )
      )
    }
    limits <- limits %||% range(values, na.rm = TRUE)
    resolved <- scales::gradient_n_pal(colors)(
      scales::rescale(
        scales::squish(values, range = limits),
        to = c(0, 1),
        from = limits
      )
    )
    resolved[is.na(resolved)] <- na_color
    resolved <- .cell_colors_to_hex(resolved)

    return(list(
      colors = colors,
      limits = limits,
      na_color = na_color,
      type = "continuous",
      resolved = resolved
    ))
  }

  character_values <- as.character(values)
  if (is.null(limits)) {
    limits <- .cell_categorical_levels(values)
  }
  if (length(limits) == 0) {
    return(list(
      colors = colors,
      limits = limits,
      na_color = na_color,
      type = "categorical",
      resolved = rep(.cell_colors_to_hex(na_color), length(values))
    ))
  }

  if (!is.null(names(colors))) {
    level_colors <- colors[limits]
  } else if (length(colors) == 1) {
    level_colors <- rep(colors, length(limits))
    names(level_colors) <- limits
  } else {
    level_colors <- grDevices::colorRampPalette(colors)(length(limits))
    names(level_colors) <- limits
  }

  resolved <- unname(level_colors[character_values])
  resolved[is.na(resolved)] <- na_color
  resolved <- .cell_colors_to_hex(resolved)

  return(list(
    colors = level_colors,
    limits = limits,
    na_color = na_color,
    type = "categorical",
    resolved = resolved
  ))
}

#' Convert opaque colors to portable hexadecimal values
#'
#' @param colors A vector of valid, fully opaque R colors.
#'
#' @return A character vector of hexadecimal colors.
#'
#' @noRd
.cell_colors_to_hex <- function(colors) {
  channels <- grDevices::col2rgb(colors)
  return(grDevices::rgb(
    red = channels[1, ],
    green = channels[2, ],
    blue = channels[3, ],
    maxColorValue = 255
  ))
}

#' Build a numeric cell plot scale
#'
#' Resolves a continuous output range or a named categorical lookup.
#'
#' @param values Mapped values or `NULL`.
#' @param specification A size or alpha specification.
#' @param output_name Name of the output-range field.
#' @param default_output Default output range.
#' @param constant Value shared by every point when no column is mapped.
#'
#' @return Resolved values and scale metadata.
#'
#' @noRd
.build_cell_numeric_scale <- function(
  values,
  specification,
  output_name,
  default_output,
  constant
) {
  if (is.null(values)) {
    scale <- list(
      output = constant,
      limits = NULL,
      values = NULL,
      resolved = constant
    )
    names(scale)[1] <- output_name
    return(scale)
  }

  mapped_values <- values
  output <- specification[[output_name]] %||% default_output
  limits <- specification$limits
  categorical <- is.character(values) || is.factor(values)
  named_output <- !is.null(names(output))

  if (categorical && named_output) {
    resolved <- unname(output[match(as.character(values), names(output))])
  } else {
    if (categorical) {
      values <- as.numeric(factor(values))
    }
    if (is.null(limits)) {
      limits <- range(values, na.rm = TRUE)
    }
    resolved <- scales::rescale(
      scales::squish(values, range = limits),
      to = output,
      from = limits
    )
  }

  scale <- list(
    output = output,
    limits = limits,
    values = mapped_values,
    resolved = resolved
  )
  names(scale)[1] <- output_name
  return(scale)
}

#' Categorical levels for a mapped column
#'
#' Factor columns keep their existing levels. Character and numeric columns use
#' sorted unique values so palette assignment does not depend on row order.
#'
#' @param values A factor, character, or numeric vector.
#'
#' @return A character vector of levels.
#'
#' @noRd
.cell_categorical_levels <- function(values) {
  if (is.factor(values)) {
    return(levels(values))
  }
  return(sort(unique(as.character(values[!is.na(values)]))))
}
