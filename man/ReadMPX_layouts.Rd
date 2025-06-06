% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_data_mpx.R
\name{ReadMPX_layouts}
\alias{ReadMPX_layouts}
\title{Load layouts from an PXL file containing MPX data}
\usage{
ReadMPX_layouts(
  filename,
  cells = NULL,
  graph_projection = c("bipartite", "Anode", "linegraph", "full"),
  verbose = TRUE
)
}
\arguments{
\item{filename}{Path to a PXL file}

\item{cells}{A character vector with component IDs. If NULL, all components are loaded.}

\item{graph_projection}{The graph projection to load. Default is 'bipartite'. If multiple
projections are present in the file, only the selected one is loaded.}

\item{verbose}{Print messages}
}
\value{
A list of lists with the layouts. At the top level, the list is split by
layout. At the second level, the list is split by component. The components are sorted in
the order they appear in the PXL file.
}
\description{
Layouts can be pre-computed with the Pixelator data processing pipeline and
are stored in a hive-styled partitioning in the PXL file. This function reads
the layouts from the PXL file and returns them as a list. Use
\code{\link{inspect_pxl_file}} to check the contents of a PXL file.
}
