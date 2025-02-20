#' Get density of points in 2 dimensions.
#' @noRd
.get2Ddensity <- function(x, y, n = 500, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg MASS} must be installed to use this function")
  }

  dens_grid <- MASS::kde2d(x, y, n = n, ...)
  x_idx <- findInterval(x, dens_grid$x)
  y_idx <- findInterval(y, dens_grid$y)
  return(dens_grid$z[cbind(x_idx, y_idx)])
}

#' @noRd
.validateDensityInputs <- function(object, marker1, marker2, facet_vars,
                                   plot_gate, gate_type = NULL,
                                   margin_density, coord_fixed) {
  # Basic object validation
  if (!inherits(object, "Seurat")) {
    cli::cli_abort("{.code object} must be a Seurat object")
  }

  # Marker validation
  if (!marker1 %in% rownames(object)) {
    cli::cli_abort("{.val {marker1}} not found in Seurat object")
  }
  if (!marker2 %in% rownames(object)) {
    cli::cli_abort("{.val {marker2}} not found in Seurat object")
  }

  # Facet variable checks
  if (!is.null(facet_vars)) {
    missing_facets <- facet_vars[!facet_vars %in% colnames(object[[]])]
    if (length(missing_facets) > 0) {
      cli::cli_abort(c(
        "Missing facet variables in Seurat object:",
        "x" = "{.val {missing_facets}}"
      ))
    }

    if (length(facet_vars) > 2) {
      cli::cli_abort("Maximum 2 facet variables supported")
    }
  }

  # Gate validation
  if (!is.null(plot_gate)) {
    if (is.null(gate_type)) {
      cli::cli_abort(
        "{.code gate_type} must be specified when using {.code plot_gate}"
      )
    }

    # Validate gate type with proper default handling
    gate_type <- match.arg(gate_type, choices = c("rectangle", "quadrant"))

    required_cols <- switch(gate_type,
                            "rectangle" = c("xmin", "xmax", "ymin", "ymax"),
                            "quadrant" = c("x", "y")
    )

    missing_cols <- setdiff(required_cols, colnames(plot_gate))
    if (length(missing_cols) > 0) {
      cli::cli_abort(c(
        "Missing required columns in {.code plot_gate}:",
        "x" = "{.val {missing_cols}}"
      ))
    }
  }

  # Compatibility checks
  if (!is.null(facet_vars) && isTRUE(margin_density)) {
    cli::cli_alert_warning(
      "Marginal density plots disabled when using faceting"
    )
    margin_density <- FALSE
  }

  if (isTRUE(coord_fixed) && isTRUE(margin_density)) {
    cli::cli_alert_warning(
      "Disabling fixed coordinates for marginal density compatibility"
    )
    coord_fixed <- FALSE
  }

  return(list(
    margin_density = margin_density,
    coord_fixed = coord_fixed
  ))
}

#' @noRd
.prepareDensityData <- function(object, marker1, marker2, facet_vars,
                                grid_n, scale_density, layer, ...) {
  plot_data <- Seurat::FetchData(
    object,
    vars = c(facet_vars, marker1, marker2),
    layer = layer
  ) %>%
    dplyr::rename(
      marker1 = !!marker1,
      marker2 = !!marker2
    )

  group_vars <- if (!is.null(facet_vars)) facet_vars else character(0)

  plot_data <- plot_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::mutate(dens = .get2Ddensity(marker1, marker2, n = grid_n, ...)) %>%
    dplyr::ungroup()

  if (scale_density) {
    plot_data <- plot_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::mutate(dens = dens / max(dens, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }

  return(plot_data)
}

#' @noRd
.createBasePlot <- function(plot_data, pt_size, alpha, colors, coord_fixed) {
  gg <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = marker1, y = marker2, color = dens)
  ) +
    ggplot2::geom_point(size = pt_size, alpha = alpha, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, color = "gray50") +
    ggplot2::geom_vline(xintercept = 0, color = "gray50") +
    ggplot2::labs(x = "Marker1", y = "Marker2") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

  if (!is.null(colors)) {
    gg <- gg + ggplot2::scale_color_gradientn(colors = colors)
  } else {
    gg <- gg + ggplot2::scale_color_viridis_c(option = "turbo")
  }

  if (coord_fixed) {
    gg <- gg + ggplot2::coord_fixed()
  }

  return(gg)
}


#' @noRd
.addGate <- function(gg, plot_gate, gate_type, facet_vars, annotation_params) {
  # Validate gate type
  gate_type <- match.arg(gate_type, c("rectangle", "quadrant"))

  # Join gate data with facet metadata
  if (!is.null(facet_vars)) {
    join_vars <- intersect(facet_vars, names(plot_gate))
    if (length(join_vars) > 0) {
      gate_data <- plot_gate %>%
        dplyr::inner_join(
          unique(gg$data[, facet_vars, drop = FALSE]),
          by = join_vars
        )
    } else {
      gate_data <- plot_gate %>%
        dplyr::cross_join(unique(gg$data[, facet_vars, drop = FALSE]))
    }
  } else {
    gate_data <- plot_gate
  }

  # Common gate validation
  if (nrow(gate_data) == 0) {
    cli::cli_abort("No matching facets found between gate data and plot data")
  }

  # Group data by facet variables if they exist
  plot_data <- gg$data
  if (!is.null(facet_vars)) {
    plot_data <- plot_data %>%
      dplyr::group_by(!!!rlang::syms(facet_vars))
  }

  # Calculate points inside gates for each facet group
  gate_data <- gate_data %>%
    dplyr::group_by(!!!rlang::syms(facet_vars)) %>%
    dplyr::group_modify(function(.x, .y) {
      # Get the current facet's data by filtering based on group keys
      current_data <- if (!is.null(facet_vars)) {
        filter_conds <- lapply(facet_vars, function(var) {
          quo <- rlang::sym(var)
          rlang::quo(!!quo == .y[[var]])
        })
        plot_data %>%
          dplyr::filter(!!!filter_conds)
      } else {
        plot_data
      }

      if (gate_type == "rectangle") {
        .x %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            n_inside = sum(
              current_data$marker1 >= xmin & current_data$marker1 <= xmax &
                current_data$marker2 >= ymin & current_data$marker2 <= ymax
            ),
            total = nrow(current_data)
          )
      } else {
        .x %>%
          tidyr::crossing(
            quadrant = c("top_left", "top_right", "bottom_left", "bottom_right")
          ) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            total = nrow(current_data),
            n_inside = sum(
              if (quadrant == "top_left") {
                current_data$marker1 < x & current_data$marker2 > y
              } else if (quadrant == "top_right") {
                current_data$marker1 >= x & current_data$marker2 > y
              } else if (quadrant == "bottom_left") {
                current_data$marker1 < x & current_data$marker2 <= y
              } else {
                current_data$marker1 >= x & current_data$marker2 <= y
              }
            )
          )
      }
    }) %>%
    dplyr::ungroup()

  # Gate-type specific processing
  if (gate_type == "rectangle") {
    # Rectangle gate validation
    req_cols <- c("xmin", "xmax", "ymin", "ymax")
    missing_cols <- setdiff(req_cols, colnames(gate_data))
    if (length(missing_cols) > 0) {
      cli::cli_abort(c(
        "Missing columns for rectangle gate:",
        "x" = "{.val {missing_cols}}"
      ))
    }

    # Calculate annotation positions
    gate_labels <- gate_data %>%
      dplyr::mutate(
        x = (xmin + xmax) / 2,
        y = ymax,
        label = sprintf("%.1f%%", 100 * (n_inside / total))
      )

    # Create gate layers
    gg <- gg +
      ggplot2::geom_rect(
        data = gate_data,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        color = "black",
        fill = NA,
        linetype = "dashed"
      )

  } else { # quadrant gate
    # Quadrant gate validation
    req_cols <- c("x", "y")
    missing_cols <- setdiff(req_cols, colnames(gate_data))
    if (length(missing_cols) > 0) {
      cli::cli_abort(c(
        "Missing columns for quadrant gate:",
        "x" = "{.val {missing_cols}}"
      ))
    }

    # Calculate plot ranges
    x_range <- ggplot2::layer_scales(gg)$x$range$range
    y_range <- ggplot2::layer_scales(gg)$y$range$range

    # Create gate labels with positions
    gate_labels <- gate_data %>%
      dplyr::mutate(
        x_label = dplyr::case_when(
          quadrant %in% c("top_left", "bottom_left") ~ x_range[1] + 0.05 * diff(x_range),
          quadrant %in% c("top_right", "bottom_right") ~ x_range[2] - 0.05 * diff(x_range)
        ),
        # Calculate consistent spacing based on plot range
        label_spacing = 0.05 * diff(y_range),
        y_label = dplyr::case_when(
          quadrant %in% c("top_left", "top_right") ~ y_range[2] - label_spacing,
          quadrant %in% c("bottom_left", "bottom_right") ~ y - label_spacing
        ),
        hjust = ifelse(quadrant %in% c("top_left", "bottom_left"), 0, 1),
        vjust = dplyr::case_when(quadrant %in% c("top_left", "top_right") ~ 0),
        label = sprintf("%.1f%%", 100 * (n_inside / total))
      )

    # Create gate layers
    gg <- gg +
      ggplot2::geom_vline(
        data = unique(gate_data[, c("x", facet_vars)]),
        ggplot2::aes(xintercept = x),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_hline(
        data = unique(gate_data[, c("y", facet_vars)]),
        ggplot2::aes(yintercept = y),
        linetype = "dashed",
        color = "black"
      )
  }

  # Add text annotations
  annotation_args <- list(
    data = gate_labels,
    inherit.aes = FALSE,
    show.legend = FALSE
  )

  # Gate-specific aesthetic mappings
  if (gate_type == "rectangle") {
    annotation_args$mapping <- ggplot2::aes(x = x, y = y, label = label)
    default_params <- list(color = "black", size = 3, vjust = 1)
  } else {
    annotation_args$mapping <- ggplot2::aes(
      x = x_label, y = y_label,
      label = label, hjust = hjust, vjust = vjust
    )
    default_params <- list(color = "black", size = 3.5)
  }

  # Merge user parameters with defaults
  final_params <- utils::modifyList(
    default_params,
    annotation_params %||% list()
  )

  # Add annotations to plot
  gg <- gg +
    do.call(ggplot2::geom_text, c(annotation_args, final_params))

  return(gg)
}


#' @noRd
.addMarginalDensity <- function(gg, plot_data) {
  x_dens <- ggplot2::ggplot(plot_data, ggplot2::aes(x = marker1)) +
    ggplot2::geom_density(fill = "grey70", color = NA) +
    ggplot2::theme_void() +
    ggplot2::scale_x_continuous(limits = range(plot_data$marker1))

  y_dens <- ggplot2::ggplot(plot_data, ggplot2::aes(x = marker2)) +
    ggplot2::geom_density(fill = "grey70", color = NA) +
    ggplot2::theme_void() +
    ggplot2::coord_flip() +
    ggplot2::scale_x_continuous(limits = range(plot_data$marker2))

  patchwork::wrap_plots(
    x_dens, patchwork::plot_spacer(), gg, y_dens,
    ncol = 2, widths = c(4, 1), heights = c(1, 4)
  )
}

#' Updated main function
#' @export
DensityScatterPlot <- function(
    object,
    marker1,
    marker2,
    facet_vars = NULL,
    plot_gate = NULL,
    gate_type = c("rectangle", "quadrant"),
    grid_n = 500,
    scale_density = TRUE,
    margin_density = TRUE,
    pt_size = 1,
    alpha = 1,
    layer = NULL,
    coord_fixed = TRUE,
    annotation_params = NULL,
    colors = NULL,
    ...) {

  # Validate inputs
  validated <- .validateDensityInputs(
    object, marker1, marker2, facet_vars,
    plot_gate, gate_type, margin_density,
    coord_fixed
  )

  # Prepare data
  plot_data <- .prepareDensityData(
    object, marker1, marker2, facet_vars,
    grid_n, scale_density, layer, ...
  )

  # Create base plot
  gg <- .createBasePlot(
    plot_data, pt_size, alpha, colors,
    validated$coord_fixed
  )

  # Add faceting
  if (!is.null(facet_vars)) {
    if (length(facet_vars) == 1) {
      gg <- gg + ggplot2::facet_grid(rows = vars(!!rlang::sym(facet_vars)))
    } else {
      gg <- gg + ggplot2::facet_grid(
        rows = vars(!!rlang::sym(facet_vars[1])),
        cols = vars(!!rlang::sym(facet_vars[2]))
      )
    }
  }

  # Add gates
  if (!is.null(plot_gate)) {
    gate_type <- match.arg(gate_type, c("rectangle", "quadrant"))
    gg <- .addGate(
      gg = gg,
      plot_gate = plot_gate,
      gate_type = gate_type,
      facet_vars = facet_vars,
      annotation_params = annotation_params
    )
  }

  # Add marginal density
  if (validated$margin_density) {
    gg <- .addMarginalDensity(gg, plot_data)
  }

return(gg)
}
