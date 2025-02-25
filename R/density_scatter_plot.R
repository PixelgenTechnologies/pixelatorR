#' Create a density scatter / pseudocolor plot.
#'
#' Create a density scatter plot, also known as a pseudocolor plot, of two markers from a Seurat object. The points in
#' the scatter plot are colored by the point density, determined by 2D kernel density estimation. Optionally, the plot
#' can be faceted by one or two variables and a gate can be plotted.
#'
#'
#' @param object A Seurat object.
#' @param marker1 Name of first marker to plot.
#' @param marker2 Name of second marker to plot.
#' @param facet_vars Optional character vector of length 1 or 2 specifying variables to facet by.
#' @param plot_gate Optional data frame containing gate parameters. For `gate_type = "rectangle"`, must contain columns
#'   'xmin', 'xmax', 'ymin', 'ymax' defining the rectangular gate boundaries. For gate_type = "quadrant", must contain
#'   columns 'x' and 'y' defining the position of the quadrant lines. The data frame can also contain the variables
#'   specified in 'facet_vars' to plot different gates in different facets.
#' @param gate_type Optional string specifying the gate type. Must be either "rectangle" or "quadrant".
#' Required if `plot_gate` is provided.
#' @param grid_n Number of grid points in each direction to compute density.
#' @param scale_density Whether to scale the density within each facet.
#' @param margin_density Whether to plot density plots in the margins.
#' @param pt_size Point size.
#' @param alpha Point transparency.
#' @param layer Optional string specifying a layer to plot.
#' @param coord_fixed Whether to use fixed coordinate ratio.
#' @param annotation_params Optional list of parameters to pass to `geom_text` for gate annotations. Common parameters
#' include color (text color), vjust (vertical justification), hjust (horizontal justification), and size (text size).
#' @param colors Optional character vector of colors to use for the color scale.
#' @param ... Additional arguments to pass to `MASS::kde2d`.
#'
#' @return A ggplot object.
#'
#'
#' @examples
#' library(pixelatorR)
#' library(Seurat)
#'
#' set.seed(123)
#'
#' # A mock-up Seurat Object
#' object <-
#'   CreateSeuratObject(counts = matrix(
#'     c(
#'       rpois(100000, 40),
#'       rpois(100000, 5)
#'     )[sample(1:200000, 200000)],
#'     nrow = 100, ncol = 2000,
#'     dimnames = list(
#'       paste0("Feature", 1:100),
#'       paste0("Cell", 1:2000)
#'     )
#'   ))
#'
#' object <-
#'   AddMetaData(object,
#'     metadata = data.frame(
#'       sample = rep(c("A", "B"), each = 1000),
#'       sample_type = rep(c("Unstimulated", "Stimulated"),
#'         each = 500, times = 2
#'       ),
#'       row.names = paste0("Cell", 1:2000)
#'     )
#'   )
#'
#' plot_gate <-
#'   data.frame(
#'     xmin = c(20, 20),
#'     xmax = c(70, 70),
#'     ymin = c(20, 20),
#'     ymax = c(60, 60),
#'     sample = c("A", "B")
#'   )
#'
#'
#' DensityScatterPlot(object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2",
#'   facet_vars = "sample",
#'   plot_gate = plot_gate,
#'   gate_type = "rectangle",
#'   layer = "counts"
#' )
#'
#'
#' # Create a Seurat object with random data
#' counts <- matrix(rpois(200, 5), nrow = 10)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Cell", 1:20)
#' object <- CreateSeuratObject(counts = counts)
#'
#' # Create a basic density scatter plot
#' DensityScatterPlot(
#'   object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2"
#' )
#'
#' # Create a density scatter plot with a quadrant gate
#' # Define quadrant gate parameters
#' quad_gate <- data.frame(
#'   x = 2, # x-coordinate for vertical line
#'   y = 3 # y-coordinate for horizontal line
#' )
#'
#' # Plot with quadrant gate and custom annotations
#' DensityScatterPlot(
#'   object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2",
#'   plot_gate = quad_gate,
#'   gate_type = "quadrant",
#'   annotation_params = list(
#'     color = "darkblue",
#'     size = 4,
#'     fontface = "bold"
#'   )
#' )
#'
#' @export
#'
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
  ...
) {
  # Validate inputs
  validated_inputs <- .validateDensityInputs(
    object,
    marker1,
    marker2,
    facet_vars,
    plot_gate,
    gate_type,
    margin_density,
    coord_fixed,
    grid_n,
    scale_density,
    pt_size,
    alpha,
    layer,
    annotation_params
  )

  # Prepare data
  plot_data <- .prepareDensityData(
    object,
    marker1,
    marker2,
    facet_vars,
    grid_n,
    scale_density,
    layer,
    ...
  )

  # Create base plot
  gg <- .createBasePlot(
    plot_data,
    pt_size,
    alpha,
    colors,
    validated_inputs$coord_fixed
  )

  # Add faceting
  if (!is.null(facet_vars)) {
    if (length(facet_vars) == 1) {
      gg <- gg + facet_grid(rows = vars(!!sym(facet_vars)))
    } else {
      gg <- gg +
        facet_grid(
          rows = vars(!!sym(facet_vars[1])),
          cols = vars(!!sym(facet_vars[2]))
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
  if (validated_inputs$margin_density) {
    gg <- .addMarginalDensity(gg, plot_data)
  }

  return(gg)
}

#' Validate inputs for density scatter plot
#'
#' @param object A Seurat object
#' @param marker1 Name of first marker
#' @param marker2 Name of second marker
#' @param facet_vars Variables to use for faceting
#' @param plot_gate Gate information for plotting
#' @param gate_type Type of gate ("rectangle" or "quadrant")
#' @param margin_density Whether to add marginal density plots
#' @param coord_fixed Whether to use fixed coordinate ratio
#' @param grid_n Number of grid points for density estimation
#' @param scale_density Whether to scale density values
#' @param pt_size Point size for plotting
#' @param alpha Point transparency
#' @param layer Layer to fetch data from
#' @param annotation_params Parameters for gate annotations
#' @return List with validated margin_density and coord_fixed values
#'
#' @noRd
#' @keywords internal
#'
.validateDensityInputs <- function(
  object,
  marker1,
  marker2,
  facet_vars,
  plot_gate,
  gate_type = NULL,
  margin_density,
  coord_fixed,
  grid_n,
  scale_density,
  pt_size,
  alpha,
  layer = NULL,
  annotation_params = NULL
) {
  # Basic object validation
  assert_class(object, "Seurat")

  # Marker validation
  assert_single_value(marker1, type = "string")
  assert_single_value(marker2, type = "string")
  assert_x_in_y(marker1, rownames(object))
  assert_x_in_y(marker2, rownames(object))

  # Numeric parameter validation
  assert_single_value(grid_n, type = "integer")
  assert_single_value(pt_size, type = "numeric")
  assert_single_value(alpha, type = "numeric")

  # Boolean parameter validation
  assert_single_value(scale_density, type = "bool")
  assert_single_value(margin_density, type = "bool")
  assert_single_value(coord_fixed, type = "bool")

  # Layer validation
  assert_single_value(layer, type = "string", allow_null = TRUE)

  # Facet variable checks
  if (!is.null(facet_vars)) {
    assert_vector(facet_vars, type = "character", n = 1, allow_null = TRUE)
    assert_max_length(facet_vars, n = 2, allow_null = TRUE)

    # Check each facet variable exists in the data
    for (facet_var in facet_vars) {
      assert_col_in_data(facet_var, object[[]], allow_null = TRUE)
    }
  }

  # Gate and annotation validation
  if (!is.null(annotation_params)) {
    assert_class(annotation_params, "list")
  }

  if (!is.null(plot_gate)) {
    assert_class(plot_gate, c("data.frame", "tbl_df"))

    if (is.null(gate_type)) {
      cli::cli_abort(
        "{.code gate_type} must be specified when {.code plot_gate} is provided"
      )
    }

    if (!gate_type[1] %in% c("rectangle", "quadrant")) {
      cli::cli_abort(c(
        "Invalid gate type provided.",
        "i" = "'gate_type' must be one of: {.code rectangle}, {.code quadrant}",
        "x" = "You provided: {.val {gate_type[1]}}"
      ))
    }

    gate_type <- match.arg(gate_type, choices = c("rectangle", "quadrant"))

    # Validate required columns based on gate type
    if (gate_type == "rectangle") {
      assert_col_in_data("xmin", plot_gate)
      assert_col_in_data("xmax", plot_gate)
      assert_col_in_data("ymin", plot_gate)
      assert_col_in_data("ymax", plot_gate)
    } else {
      assert_col_in_data("x", plot_gate)
      assert_col_in_data("y", plot_gate)
    }

    # Check for invalid columns in plot_gate
    allowed_cols <- switch(
      gate_type,
      "rectangle" = c("xmin", "xmax", "ymin", "ymax"),
      "quadrant" = c("x", "y")
    )
    check <- !names(plot_gate) %in% c(allowed_cols, facet_vars)
    if (sum(check) > 0) {
      cli::cli_abort(
        c(
          "i" = "{.code plot_gate} can only contain gate-specific columns and faceting variables",
          "x" = "The following columns are not allowed: ",
          " " = "{.val {names(plot_gate)[check]}}"
        )
      )
    }
  }

  # Compatibility checks
  if (!is.null(facet_vars) && isTRUE(margin_density)) {
    cli::cli_alert_warning(
      "Marginal density ({.code margin_density = TRUE}) is not supported with faceting.
      Setting {.code margin_density = FALSE}"
    )
    margin_density <- FALSE
  }

  if (isTRUE(coord_fixed) && isTRUE(margin_density)) {
    cli::cli_alert_warning(
      "Setting {.code coord_fixed = FALSE} for compatibility with {.code margin_density = TRUE}"
    )
    coord_fixed <- FALSE
  }

  return(list(
    margin_density = margin_density,
    coord_fixed = coord_fixed
  ))
}

#' Get density of points in 2 dimensions.
#'
#' Calculate the kernel density of points in 2 dimensions.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @param ... Additional arguments to pass to 'MASS::kde2d'.
#' @return The density within each square.
#'
#' @noRd
#' @keywords internal
#'
.get2Ddensity <- function(x, y, n = 500, ...) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg MASS} must be installed to use this function")
  }

  dens_grid <- MASS::kde2d(x, y, n = n, ...)
  x_idx <- findInterval(x, dens_grid$x)
  y_idx <- findInterval(y, dens_grid$y)
  return(dens_grid$z[cbind(x_idx, y_idx)])
}

#' Prepare data for density scatter plot
#'
#' @param object A Seurat object
#' @param marker1 Name of first marker
#' @param marker2 Name of second marker
#' @param facet_vars Variables to use for faceting
#' @param grid_n Size of density estimation grid
#' @param scale_density Whether to scale density values
#' @param layer Which layer to fetch data from
#' @param ... Additional arguments passed to .get2Ddensity
#' @return A data frame with density values
#'
#' @noRd
#' @keywords internal
.prepareDensityData <- function(
  object,
  marker1,
  marker2,
  facet_vars,
  grid_n,
  scale_density,
  layer,
  ...
) {
  plot_data <- FetchData(
    object,
    vars = c(facet_vars, marker1, marker2),
    layer = layer
  ) %>%
    rename(
      marker1 = !!marker1,
      marker2 = !!marker2
    )

  group_vars <- if (!is.null(facet_vars)) facet_vars else character(0)

  plot_data <- plot_data %>%
    group_by(across(all_of(group_vars))) %>%
    mutate(dens = .get2Ddensity(marker1, marker2, n = grid_n, ...)) %>%
    ungroup()

  if (scale_density) {
    plot_data <- plot_data %>%
      group_by(across(all_of(group_vars))) %>%
      mutate(dens = dens / max(dens, na.rm = TRUE)) %>%
      ungroup()
  }

  return(plot_data)
}

#' Create base ggplot for density scatter plot
#'
#' @param plot_data Data frame with plot data
#' @param pt_size Point size for scatter plot
#' @param alpha Alpha transparency for points
#' @param colors Vector of colors for density gradient
#' @param coord_fixed Whether to use fixed coordinate ratio
#' @return A ggplot object
#'
#' @noRd
#' @keywords internal
#'
.createBasePlot <- function(plot_data, pt_size, alpha, colors, coord_fixed) {
  gg <- ggplot(
    plot_data,
    aes(x = marker1, y = marker2, color = dens)
  ) +
    geom_point(size = pt_size, alpha = alpha, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_vline(xintercept = 0, color = "gray50") +
    labs(x = "Marker1", y = "Marker2") +
    theme_bw() +
    theme(panel.grid = element_blank())

  if (!is.null(colors)) {
    gg <- gg + scale_color_gradientn(colors = colors)
  } else {
    gg <- gg + scale_color_viridis_c(option = "turbo")
  }

  if (coord_fixed) {
    gg <- gg + coord_fixed()
  }

  return(gg)
}

#' Add gate annotation to density scatter plot
#'
#' @param gg Base ggplot object
#' @param plot_gate Gate information for plotting
#' @param gate_type Type of gate ("rectangle" or "quadrant")
#' @param facet_vars Variables used for faceting
#' @param annotation_params Parameters for gate annotation
#' @return A ggplot object with gate annotation
#'
#' @noRd
#' @keywords internal
#'
.addGate <- function(gg, plot_gate, gate_type, facet_vars, annotation_params) {
  # Validate gate type
  gate_type <- match.arg(gate_type, c("rectangle", "quadrant"))

  # Join gate data with facet metadata
  if (!is.null(facet_vars)) {
    join_vars <- intersect(facet_vars, names(plot_gate))
    if (length(join_vars) > 0) {
      gate_data <- plot_gate %>%
        inner_join(
          unique(gg$data[, facet_vars, drop = FALSE]),
          by = join_vars
        )
    } else {
      gate_data <- plot_gate %>%
        cross_join(unique(gg$data[, facet_vars, drop = FALSE]))
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
      group_by(!!!syms(facet_vars))
  }

  # Calculate points inside gates for each facet group
  gate_data <- gate_data %>%
    group_by(!!!syms(facet_vars)) %>%
    group_modify(function(.x, .y) {
      # Get the current facet's data by filtering based on group keys
      current_data <- if (!is.null(facet_vars)) {
        filter_conds <- lapply(facet_vars, function(var) {
          facet_var <- sym(var)
          quo(!!facet_var == .y[[var]])
        })
        plot_data %>%
          filter(!!!filter_conds)
      } else {
        plot_data
      }

      if (gate_type == "rectangle") {
        .x %>%
          rowwise() %>%
          mutate(
            n_inside = sum(
              current_data$marker1 >= xmin &
                current_data$marker1 <= xmax &
                current_data$marker2 >= ymin &
                current_data$marker2 <= ymax
            ),
            total = nrow(current_data)
          )
      } else {
        .x %>%
          crossing(
            quadrant = c("top_left", "top_right", "bottom_left", "bottom_right")
          ) %>%
          rowwise() %>%
          mutate(
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
    ungroup()

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
      mutate(
        x = (xmin + xmax) / 2,
        y = ymax,
        label = sprintf("%.1f%%", 100 * (n_inside / total))
      )

    # Create gate layers
    gg <- gg +
      geom_rect(
        data = gate_data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        color = "black",
        fill = NA,
        linetype = "dashed"
      )
  } else {
    # quadrant gate
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
    x_range <- layer_scales(gg)$x$range$range
    y_range <- layer_scales(gg)$y$range$range

    # Create gate labels with positions
    gate_labels <- gate_data %>%
      mutate(
        x_label = case_when(
          quadrant %in% c("top_left", "bottom_left") ~
            x_range[1] + 0.05 * diff(x_range),
          quadrant %in% c("top_right", "bottom_right") ~
            x_range[2] - 0.05 * diff(x_range)
        ),
        # Calculate consistent spacing based on plot range
        label_spacing = 0.05 * diff(y_range),
        y_label = case_when(
          quadrant %in% c("top_left", "top_right") ~ y_range[2] - label_spacing,
          quadrant %in% c("bottom_left", "bottom_right") ~ y - label_spacing
        ),
        hjust = ifelse(quadrant %in% c("top_left", "bottom_left"), 0, 1),
        vjust = case_when(quadrant %in% c("top_left", "top_right") ~ 0),
        label = sprintf("%.1f%%", 100 * (n_inside / total))
      )

    # Create gate layers
    gg <- gg +
      geom_vline(
        data = unique(gate_data[, c("x", facet_vars)]),
        aes(xintercept = x),
        linetype = "dashed",
        color = "black"
      ) +
      geom_hline(
        data = unique(gate_data[, c("y", facet_vars)]),
        aes(yintercept = y),
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
    annotation_args$mapping <- aes(x = x, y = y, label = label)
    default_params <- list(color = "black", size = 3, vjust = 1)
  } else {
    annotation_args$mapping <- aes(
      x = x_label,
      y = y_label,
      label = label,
      hjust = hjust,
      vjust = vjust
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
    do.call(geom_text, c(annotation_args, final_params))

  return(gg)
}

#' Add marginal density plots
#'
#' @param gg Base ggplot object
#' @param plot_data Data frame with plot data
#' @return A ggplot object with marginal density plots
#'
#' @noRd
#' @keywords internal
#'
.addMarginalDensity <- function(gg, plot_data) {
  x_dens <- ggplot(plot_data, aes(x = marker1)) +
    geom_density(fill = "grey70", color = NA) +
    theme_void() +
    scale_x_continuous(limits = range(plot_data$marker1))

  y_dens <- ggplot(plot_data, aes(x = marker2)) +
    geom_density(fill = "grey70", color = NA) +
    theme_void() +
    coord_flip() +
    scale_x_continuous(limits = range(plot_data$marker2))

  wrap_plots(
    x_dens,
    plot_spacer(),
    gg,
    y_dens,
    ncol = 2,
    widths = c(4, 1),
    heights = c(1, 4)
  )
}

#' Create a density scatter / pseudocolor plot.
#'
#' Create a density scatter plot, also known as a pseudocolor plot, of two markers from a Seurat object. The points in
#' the scatter plot are colored by the point density, determined by 2D kernel density estimation. Optionally, the plot
#' can be faceted by one or two variables and a gate can be plotted.
#'
#'
#' @param object A Seurat object.
#' @param marker1 Name of first marker to plot.
#' @param marker2 Name of second marker to plot.
#' @param facet_vars Optional character vector of length 1 or 2 specifying variables to facet by.
#' @param plot_gate Optional data frame containing gate parameters. If provided, `gate_type` must also be specified.
#' @param gate_type Optional string specifying the gate type. Must be either "rectangle" or "quadrant".
#' Required if `plot_gate` is provided.
#' @param grid_n Number of grid points in each direction to compute density.
#' @param scale_density Whether to scale the density within each facet.
#' @param margin_density Whether to plot density plots in the margins.
#' @param pt_size Point size.
#' @param alpha Point transparency.
#' @param layer Optional string specifying a layer to plot.
#' @param coord_fixed Whether to use fixed coordinate ratio.
#' @param annotation_params Optional list of parameters to pass to `geom_text` for gate annotations. Common parameters
#' include color (text color), vjust (vertical justification), hjust (horizontal justification), and size (text size).
#' @param colors Optional character vector of colors to use for the color scale.
#' @param ... Additional arguments to pass to `MASS::kde2d`.
#'
#' @return A ggplot object.
#'
#'
#' @examples
#'
#' library(pixelatorR)
#' library(Seurat)
#'
#' set.seed(123)
#'
#' # A mock-up Seurat Object
#' object <-
#'   CreateSeuratObject(counts = matrix(
#'     c(
#'       rpois(100000, 40),
#'       rpois(100000, 5)
#'     )[sample(1:200000, 200000)],
#'     nrow = 100, ncol = 2000,
#'     dimnames = list(
#'       paste0("Feature", 1:100),
#'       paste0("Cell", 1:2000)
#'     )
#'   ))
#'
#' object <-
#'   AddMetaData(object,
#'     metadata = data.frame(
#'       sample = rep(c("A", "B"), each = 1000),
#'       sample_type = rep(c("Unstimulated", "Stimulated"),
#'         each = 500, times = 2
#'       ),
#'       row.names = paste0("Cell", 1:2000)
#'     )
#'   )
#'
#' plot_gate <-
#'   data.frame(
#'     xmin = c(20, 20),
#'     xmax = c(70, 70),
#'     ymin = c(20, 20),
#'     ymax = c(60, 60),
#'     sample = c("A", "B")
#'   )
#'
#'
#' DensityScatterPlot(object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2",
#'   facet_vars = "sample",
#'   plot_gate = plot_gate,
#'   gate_type = "rectangle",
#'   layer = "counts"
#' )
#'
#'
#' # Create a Seurat object with random data
#' counts <- matrix(rpois(200, 5), nrow = 10)
#' rownames(counts) <- paste0("Feature", 1:10)
#' colnames(counts) <- paste0("Cell", 1:20)
#' object <- CreateSeuratObject(counts = counts)
#'
#' # Create a basic density scatter plot
#' DensityScatterPlot(
#'   object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2"
#' )
#'
#' # Create a density scatter plot with a quadrant gate
#' # Define quadrant gate parameters
#' quad_gate <- data.frame(
#'   x = 2, # x-coordinate for vertical line
#'   y = 3 # y-coordinate for horizontal line
#' )
#'
#' # Plot with quadrant gate and custom annotations
#' DensityScatterPlot(
#'   object,
#'   marker1 = "Feature1",
#'   marker2 = "Feature2",
#'   plot_gate = quad_gate,
#'   gate_type = "quadrant",
#'   annotation_params = list(
#'     color = "darkblue",
#'     size = 4,
#'     fontface = "bold"
#'   )
#' )
#'
#' @export
#'
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
  ...
) {
  # Validate inputs
  validated_inputs <- .validateDensityInputs(
    object,
    marker1,
    marker2,
    facet_vars,
    plot_gate,
    gate_type,
    margin_density,
    coord_fixed,
    grid_n,
    scale_density,
    pt_size,
    alpha,
    layer,
    annotation_params
  )

  # Prepare data
  plot_data <- .prepareDensityData(
    object,
    marker1,
    marker2,
    facet_vars,
    grid_n,
    scale_density,
    layer,
    ...
  )

  # Create base plot
  gg <- .createBasePlot(
    plot_data,
    pt_size,
    alpha,
    colors,
    validated_inputs$coord_fixed
  )

  # Add faceting
  if (!is.null(facet_vars)) {
    if (length(facet_vars) == 1) {
      gg <- gg + facet_grid(rows = vars(!!sym(facet_vars)))
    } else {
      gg <- gg +
        facet_grid(
          rows = vars(!!sym(facet_vars[1])),
          cols = vars(!!sym(facet_vars[2]))
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
  if (validated_inputs$margin_density) {
    gg <- .addMarginalDensity(gg, plot_data)
  }

  return(gg)
}
