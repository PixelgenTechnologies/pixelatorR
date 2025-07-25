% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequencing_saturation.R
\name{sequencing_saturation}
\alias{sequencing_saturation}
\title{Calculate Sequencing Saturation}
\usage{
sequencing_saturation(graph_elements, graph_reads)
}
\arguments{
\item{graph_elements}{The number of graph elements (e.g., nodes, edges).}

\item{graph_reads}{The number of reads or supporting elements.}
}
\value{
The sequencing saturation of the graph expressed as a percentage.
}
\description{
This function calculates the sequencing saturation of a sample or a graph component,
which can be applied to unique reads, nodes, or edges.
}
\details{
The sequencing saturation is calculated using the formula:
\deqn{S = 100 \times \left(1 - \frac{E}{R}\right)}
where:
\itemize{
\item \eqn{S} is the saturation (as a percentage),
\item \eqn{E} is the number of unique reads or graph elements (e.g., nodes, edges), and
\item \eqn{R} is the number of reads or supporting elements in the graph.
}

This function can be used for calculating the saturation of:
\itemize{
\item Reads: Number of deduplicated reads vs. reads
\item Nodes: Number of nodes vs. reads
\item Edges: Number of edges vs. reads
\item Other graph elements: Adjust the input accordingly
}
}
\examples{

# For a graph with 300 unique reads, 100 nodes, and 200 edges,
# sequenced at 400 total reads

# Read sequencing saturation
sequencing_saturation(300, 400)

# Node sequencing saturation
sequencing_saturation(100, 400)

# Edge sequencing saturation
sequencing_saturation(200, 400)

}
