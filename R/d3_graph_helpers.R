#' Add node colors to a \code{CellGraph}
#'
#' Adds node colors based on the expression of a single or multiple
#' markers. If multiple markers are provided, their values are either
#' multiplied or summed depending on the \code{mode}.
#'
#' @param cg A \code{\link{CellGraph}} object
#' @param markers A character vector with valid marker IDs
#' @param smooth_counts Applies neighborhood smoothing of marker counts.
#' Each node marker counts is combined with the sum of marker counts in
#' a neighborhood of length 1.
#' @param palette A valid color palette from \code{RColorBrewer} or \code{viridis}
#' or a vector of colors.
#' @param rev_pal Should the palette be reversed?
#' @param mode Only used if more than 1 marker is provided. With \code{mode = "product"},
#' the marker counts are multiplied and with \code{mode = "sum"}, the marker counts are
#' summed. This step is executed after smoothing.
#' @param normalize Apply normalization to the marker counts after smoothing and
#' the multiplication/summation step. The values are first divided by the node degree
#' followed by log-transformation with \code{log1p}.
#' @param trim_quantiles A numeric vector of length 2 specifying quantiles to trim.
#' This can be useful to reduce the influence of outliers. The specified quantiles
#' are calculated from the marker counts values prior to normalization.
#' (see section below for details)
#' @param nNodes Number of nodes to keep in the graph
#'
#' @section Trim quantile:
#' Sometimes, a few outliers can dominate which makes the colors difficult to interpret.
#' With \code{trim_quantiles}, we can remove such outliers. By setting the upper bound to
#' 0.99, we are calculating the 99th quantile of the marker count values and any values
#' above this threshold are replaced with the value for the 99th quantile. The same
#' thing applies to the lower bound but the other way around.
#'
#' @export
#'
color_by_marker <- function(
  cg,
  markers,
  smooth_counts = TRUE,
  palette = "viridis",
  rev_pal = FALSE,
  mode = c("product", "sum"),
  normalize = TRUE,
  trim_quantiles = c(0, 1),
  nNodes = NULL
) {
  # Set global variables to NULL (required by shinytest2)
  . <- val <- NULL

  # Require scales library
  expect_scales()

  # Validate input parameters
  assert_class(cg, "CellGraph")
  assert_vector(markers, type = "character", n = 1)
  assert_vector(palette, type = "character", n = 1)
  assert_single_value(smooth_counts, type = "bool")
  assert_single_value(normalize, type = "bool")
  assert_class(trim_quantiles, "numeric")
  assert_length(trim_quantiles, 2)
  assert_within_limits(trim_quantiles, c(0, 1))
  if (trim_quantiles[2] <= trim_quantiles[1]) {
    cli::cli_abort(
      c("x" = "'trim_quantiles[2]' must be larger than 'trim_quantiles[1]'")
    )
  }
  if (!is.null(nNodes)) {
    assert_class(nNodes, "numeric")
    assert_within_limits(nNodes, c(1, length(cg@cellgraph)))
  }

  mode <- match.arg(mode, choices = c("product", "sum"))

  # Fetch node counts and make sure its non-empty
  counts <- slot(cg, name = "counts")
  if (is.null(counts)) {
    cli::cli_abort(
      c("x" = "counts are missing from the 'CellGraph' object")
    )
  }

  # Validate markers
  assert_x_in_y(markers, colnames(counts))

  # Fetch tbl_graph
  g <- slot(cg, name = "cellgraph")

  # Smooth counts if smooth_counts=TRUE
  if (smooth_counts) {
    adjMat <- g %>% igraph::as_adjacency_matrix()
    counts <- adjMat %*% counts + counts
  }

  g <- g %>%
    {
      # Aggregate marker counts if more than 1 marker is provided
      if (length(markers) > 1) {
        if (mode == "product") {
          mutate(., val = apply(counts[, markers], 1, prod))
        } else if (mode == "sum") {
          mutate(., val = rowSums(counts[, markers]))
        }
      } else {
        mutate(., val = counts[, markers, drop = TRUE]) # Include marker counts in node table
      }
    } %>%
    {
      # Apply trimming if non-default trim_quantiles are provided
      if ((trim_quantiles[1] > 0) || (trim_quantiles[2] < 1)) {
        mutate(., val = .trim_quantiles(val,
          top_q = trim_quantiles[2],
          bottom_q = trim_quantiles[1]
        ))
      } else {
        .
      }
    } %>%
    {
      # Normalize values if normalize=TRUE
      if (normalize) {
        mutate(., val = log1p(val / igraph::degree(g)))
      } else {
        .
      }
    }

  # Sample graph is desired
  if (!is.null(nNodes)) {
    nNodes <- round(nNodes)
    set.seed(123)

    nodes_keep <- sample(c(rep(TRUE, nNodes), rep(FALSE, length(cg@cellgraph) - nNodes)))
    g <- g %>%
      filter(nodes_keep)

    cg@counts <- cg@counts[nodes_keep, ]
    cg@layout <- lapply(cg@layout, function(x) {
      x[nodes_keep, ]
    })
  }

  g <- g %>%
    # Compute node colors using provided palette
    mutate(color = scales::col_numeric(domain = c(0, max(val)), palette = palette, reverse = rev_pal)(val))

  slot(cg, name = "cellgraph") <- g

  return(cg)
}

#' Trim quantiles
#'
#' @param x A numeric vector
#' @param bottom_q,top_q Quantiles used for trimming
#'
#' @noRd
#'
.trim_quantiles <- function(
  x,
  bottom_q = 0,
  top_q = 0.99
) {
  # Validate x
  assert_non_empty_object(x, classes = "numeric")

  # Calculate quantiles
  low_thr <- quantile(x, probs = bottom_q)
  high_thr <- quantile(x, probs = top_q)

  # Trim data
  x[x < low_thr] <- low_thr
  x[x > high_thr] <- high_thr

  return(x)
}


#' Convert a \code{tbl_graph} to json format
#'
#' @param data A \code{tbl_graph} object
#'
#' @noRd
#'
.convert_tbl_graph_to_json <- function(
  data
) {
  # Set global variables to NULL (required by shinytest2)
  from <- to <- NULL

  # Validate data
  assert_class(data, "tbl_graph")

  # Fetch node table and add IDs
  nodes <- data %N>%
    as_tibble() %>%
    mutate(id = seq_len(n())) %>%
    mutate_if(is.factor, as.character)

  # Feth edge table and rename from/to
  links <- data %E>%
    as_tibble() %>%
    rename(source = from, target = to) %>%
    mutate_if(is.factor, as.character)

  # Assemble nodes/edges in a list and convert to json
  graph <- list(nodes = nodes, links = links)
  graph_json <- jsonlite::toJSON(graph, auto_unbox = TRUE, dataframe = "rows", pretty = TRUE)

  return(graph_json)
}

#' Add layout coordinates to nodes in a \code{tbl_graph} object
#'
#' @param cg A \code{tbl_graph} object
#' @param layout_coordinates A \code{tbl_df} with layout coordinates
#' @param scale rescale coordinates
#' @param keep_aspect_ratio Should the aspect ratio be kept?
#'
#' @noRd
#'
.add_coordinates_to_tbl_graph <- function(
  cg,
  layout_coordinates,
  scale = TRUE,
  keep_aspect_ratio = TRUE
) {
  # Validate input
  assert_class(cg, "CellGraph")
  assert_class(layout_coordinates, "data.frame")
  assert_single_value(scale, type = "bool")
  assert_single_value(keep_aspect_ratio, type = "bool")
  if (!all(c("x", "y") %in% colnames(layout_coordinates))) {
    cli::cli_abort(
      c("x" = "Columns {.str x} and {.str y} must be present in {.var layout_coordinates}")
    )
  } else {
    if ("z" %in% colnames(layout_coordinates)) {
      classes <- sapply(layout_coordinates[, c("x", "y", "z")], function(x) inherits(x, what = "numeric"))
      if (!all(classes)) {
        cli::cli_abort(
          c("x" = "Columns {.str x}, {.str y} and {.str z} must be numeric")
        )
      }
    } else {
      classes <- sapply(layout_coordinates[, c("x", "y")], function(x) inherits(x, what = "numeric"))
      if (!all(classes)) {
        cli::cli_abort(
          c("x" = "Columns {.str x} and {.str y} must be numeric")
        )
      }
    }
  }
  assert_singles_match(nrow(layout_coordinates), length(cg@cellgraph))

  # Add layout_coordinates to data@cellgraph
  cg@cellgraph <- cg@cellgraph %N>%
    mutate(
      x = layout_coordinates$x,
      y = layout_coordinates$y
    ) %>%
    {
      if ("z" %in% colnames(layout_coordinates)) {
        mutate(., z = layout_coordinates$z)
      } else {
        .
      }
    }

  # Rescale
  # nolint start
  if (scale) {
    if (keep_aspect_ratio) {
      max_val <- cg@cellgraph %>%
        as_tibble() %>%
        select(any_of(c("x", "y", "z"))) %>%
        max()
      cg@cellgraph <- cg@cellgraph %>%
        mutate(across(any_of(c("x", "y", "z")), ~ .x / max_val))
    } else {
      cg@cellgraph <- cg@cellgraph %>%
        mutate(across(any_of(c("x", "y", "z")), ~ .x / max(.x)))
    }
  }
  # nolint end

  return(cg)
}
