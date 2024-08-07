% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_layout.R
\name{layout coordinates utils}
\alias{layout coordinates utils}
\alias{center_layout_coordinates}
\alias{normalize_layout_coordinates}
\alias{project_layout_coordinates_on_unit_sphere}
\title{Layout Coordinates Utility Functions}
\usage{
center_layout_coordinates(layout)

normalize_layout_coordinates(layout)

project_layout_coordinates_on_unit_sphere(layout)
}
\arguments{
\item{layout}{A matrix-like object with layout coordinates}
}
\value{
A \code{tbl_df} object with adjusted layout coordinates
}
\description{
Utility function used to manipulate layout coordinates.
These functions always returns a \code{tbl_df} object.
}
\section{Center Layout Coordinates}{

Centers each axis of the layout coordinates around their means.
}

\section{Normalize Layout Coordinates}{

Centers each axis of the layout coordinates around the mean and
adjusts each point coordinate such that their median length
(euclidean norm) is 1.
}

\section{Project Layout Coordinates on a Unit Sphere}{

Centers each axis of the layout coordinates around the mean and
adjusts each coordinate such that their lengths (euclidean norm)
are 1. This function only accepts layouts with 3 dimensions.
}

\examples{
library(dplyr)
library(tibble)

# Generate random points that are offset to (100, 100, 100)
xyz <- matrix(rnorm(600, mean = 100, sd = 20),
  ncol = 3,
  dimnames = list(NULL, c("x", "y", "z"))
) \%>\%
  as_tibble()

# Visualize random points
plotly::plot_ly(
  data = xyz, x = ~x, y = ~y, z = ~z,
  type = "scatter3d", mode = "markers"
)

# Center points at (0, 0, 0)
xyz_centered <- center_layout_coordinates(xyz)
apply(xyz_centered, 2, mean)
plotly::plot_ly(
  data = xyz_centered, x = ~x, y = ~y,
  z = ~z, type = "scatter3d", mode = "markers"
)

# Normalize points to have a median radius of 1
xyz_normalized <- normalize_layout_coordinates(xyz_centered)
radii <- sqrt(rowSums(xyz_normalized^2))
median(radii)
plotly::plot_ly(
  data = xyz_normalized, x = ~x, y = ~y,
  z = ~z, type = "scatter3d", mode = "markers"
)

# Project points on unit sphere
xyz_projected <- project_layout_coordinates_on_unit_sphere(xyz_normalized)
radii <- sqrt(rowSums(xyz_projected^2))
all(near(radii, y = rep(1, length(radii)), tol = 1e-12))
plotly::plot_ly(
  data = xyz_projected, x = ~x, y = ~y,
  z = ~z, type = "scatter3d", mode = "markers"
)

}
