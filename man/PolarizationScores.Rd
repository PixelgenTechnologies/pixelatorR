% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/CellGraphAssay.R, R/objects.R
\name{PolarizationScores}
\alias{PolarizationScores}
\alias{PolarizationScores<-}
\alias{PolarizationScores.MPXAssay}
\alias{PolarizationScores<-.MPXAssay}
\alias{PolarizationScores.Seurat}
\alias{PolarizationScores<-.Seurat}
\title{PolarizationScores}
\usage{
PolarizationScores(object, ...)

PolarizationScores(object, ...) <- value

\method{PolarizationScores}{MPXAssay}(object, add_marker_counts = FALSE, ...)

\method{PolarizationScores}{MPXAssay}(object, ...) <- value

\method{PolarizationScores}{Seurat}(
  object,
  assay = NULL,
  meta_data_columns = NULL,
  add_marker_counts = FALSE,
  ...
)

\method{PolarizationScores}{Seurat}(object, assay = NULL, ...) <- value
}
\arguments{
\item{object}{An object with polarization scores}

\item{...}{Not implemented}

\item{value}{A \code{tbl_df} with polarization scores}

\item{add_marker_counts}{A logical value indicating whether to add marker
counts to the polarization score table.}

\item{assay}{Name of a \code{CellGraphAssay}}

\item{meta_data_columns}{A character vector with meta.data column names.
This option can be useful to join meta.data columns with the polarization
score table.}
}
\value{
\code{PolarizationScores}: Polarization scores

\code{PolarizationScores<-}: An object with polarization scores
updated
}
\description{
Get and set polarization scores for a \code{\link{CellGraphAssay}},
\code{\link{CellGraphAssay5}} or a \code{Seurat} object
}
\seealso{
Other spatial metrics: 
\code{\link{ColocalizationScores}()},
\code{\link{Edgelists}()},
\code{\link{ProximityScores}()}
}
\concept{spatial metrics}
