#' @include generics.R
NULL

#' Convert objects to a \code{\link{CellGraphAssay5}}
#'
#' @param x An object to convert to class \code{\link{CellGraphAssay5}}
#' @param ... Arguments passed to other methods
#' @rdname as.CellGraphAssay5
#' @export as.CellGraphAssay5
as.CellGraphAssay5 <- function (
  x,
  ...
) {
  UseMethod(generic = "as.CellGraphAssay5", object = x)
}

#' Convert objects to a \code{\link{CellGraphAssay}}
#'
#' @param x An object to convert to class \code{\link{CellGraphAssay}}
#' @param ... Arguments passed to other methods
#' @rdname as.CellGraphAssay
#' @export as.CellGraphAssay
as.CellGraphAssay <- function (
  x,
  ...
) {
  UseMethod(generic = "as.CellGraphAssay", object = x)
}

#' CellGraphs
#'
#' Get and set \code{\link{CellGraph}} lists for different objects.
#'
#' @param object A \code{Seurat}, \code{CellGraphAssay} or \code{CellGraphAssay} object
#' @param ... Arguments passed to other methods
#'
#' @return Returns a list of \code{\link{CellGraph}} objects. If there are
#' no \code{\link{CellGraph}} objects present, returns an empty named list.
#'
#' @rdname CellGraphs
#' @seealso [PolarizationScores()] and [ColocalizationScores()] for getting/setting spatial metrics
#'
#' @export CellGraphs
#'
CellGraphs <- function (
  object,
  ...
) {
  UseMethod(generic = "CellGraphs", object = object)
}

#' @param value A list of \code{\link{CellGraph}} objects
#'
#' @rdname CellGraphs
#' @export CellGraphs<-
#'
"CellGraphs<-" <- function (
  object,
  ...,
  value
) {
  UseMethod(generic = 'CellGraphs<-', object = object)
}

#' Load CellGraphs
#'
#' Loads a list of \code{\link{CellGraph}} objects. One can specify cells to load by
#' their names using the \code{cells} parameter.
#'
#' Graphs can be loaded as one of 'bipartite', 'Anode' or 'linegraph'. See details about each of
#' these graph representations in the sections below.
#'
#' @section Bi-partite graph:
#' In the bi-partite graph, edges can only go from a upia to a upib. The bi-partite graph
#' is first collapsed from a multigraph to a simple graph by aggregating marker counts
#' for each upia/upib combination. For visualization of marker counts on the graph, it is
#' often convenient to project values on the nodes; however, the marker counts are not available
#' in the nodes in the bi-partite graph. To circumvent this issue, node counts are calculated
#' by aggregating its edge counts. This means that the total marker count will be inflated.
#'
#' @section A node projected-graph:
#' In the A node projected-graph, pairs of upias that share a upib are connected with an edge.
#' The (upia) node counts are calculated by aggregating counts across all edges associated
#' with the A nodes. The number of nodes in the A node projected-graph is substantially
#' lower than the bi-partite graph.
#'
#' @section Line graph:
#' Starting with a bi-partite graph, a node is placed on each edge. Then, an edge is drawn
#' between each pair of adjacent nodes. The number of nodes in the linegraph is equivalent
#' to the number of edges in the bi-partite graph. However, the number of edges is substantially
#' larger. One major advantage with the linegraph is that the node counts represent the actual
#' raw data, which is not the case for the bi-partite graph and the A node projected-graph. On the
#' down side, linegraphs are much larger and tends to slow down layout computations and
#' visualizations.
#'
#' @param object A \code{Seurat} object or an \code{\link{CellGraphAssay}} object
#' @param ... Parameters passed to other methods
#'
#' @rdname LoadCellGraphs
#'
#' @return An object with a list of \code{\link{CellGraph}} objects
#'
#' @export
#'
LoadCellGraphs <- function (
  object,
  ...
) {
  UseMethod(generic = 'LoadCellGraphs', object = object)
}

#' Remove CellGraphs
#'
#' Clears the \code{\link{CellGraph}} objects from the "cellgraphs" slot.
#'
#' @param object A \code{Seurat} object or an \code{\link{CellGraphAssay}} object
#' @param ... Parameters passed to other methods
#'
#' @rdname RemoveCellGraphs
#'
#' @return An object without CellGraphs
#'
#' @export
#'
RemoveCellGraphs <- function (
  object,
  ...
) {
  UseMethod(generic = 'RemoveCellGraphs', object = object)
}

#' Edge Rank Plot
#'
#' Plots the number of edges/molecules per component against the component rank
#'
#' @param object A \code{data.frame}-like object or a \code{Seurat} object
#' @param ... Parameters passed to other methods
#'
#' @rdname EdgeRankPlot
#'
#' @return A \code{ggplot} object
#'
#' @export EdgeRankPlot
#'
EdgeRankPlot <- function (
  object,
  ...
) {
  UseMethod(generic = "EdgeRankPlot", object = object)
}

#' Plot cell counts per group
#'
#' QC plot function used to get a quick overview of called cells in an
#' MPX data set.
#'
#' @concept plots
#' @family QC-plots
#'
#' @param object A \code{data.frame}-like object or a \code{Seurat} object
#' @param group_by A column in the object representing a 'character' or 'factor'
#' to group data by
#' @param color_by A column in the object representing a 'character' or 'factor'
#' to color data by
#' @param show_count Place the count on top of the bar or next to the bar if
#' \code{flip_axes = TRUE}
#' @param flip_axes Flip the plot layout
#' @param ... Not yet implemented
#'
#' @rdname CellCountPlot
#'
#' @return A \code{ggplot} object
#'
#' @export CellCountPlot
#'
CellCountPlot <- function (
  object,
  ...
) {
  UseMethod(generic = "CellCountPlot", object = object)
}

#' Plot UMIs per UPIa for quality control
#'
#' @concept plots
#' @family QC-plots
#'
#' @param object A \code{data.frame}-like object or a \code{Seurat} object with
#' \code{umi_per_upia}, \code{tau} and \code{tau_type} values
#' @param group_by A column in the object representing a 'character' or 'factor'
#' to group data by
#' @param ... Not yet implemented
#'
#' @rdname TauPlot
#'
#' @return A \code{ggplot} object
#'
#' @export TauPlot
#'
TauPlot <- function (
  object,
  ...
) {
  UseMethod(generic = "TauPlot", object = object)
}


#' Compute a graph layout
#'
#' @param object An object
#' @param ... Additional parameters passed to other methods
#'
#' @rdname ComputeLayout
#'
#' @return An object containing a graph layout
#'
#' @export
#'
ComputeLayout <- function (
  object,
  ...
) {
  UseMethod(generic = "ComputeLayout", object = object)
}

#' Keep largest component
#'
#' Finds connected components of a graph and returns the largest component
#'
#' @param object An object
#' @param ... Parameters passed to other methods
#'
#' @rdname KeepLargestComponent
#'
#' @export
#'
KeepLargestComponent <- function (
  object,
  ...
) {
  UseMethod(generic = "KeepLargestComponent", object = object)
}

#' A-node projection
#'
#' A function to create an A node projection graph from an edgelist
#' representing a bipartite graph.
#'
#' @param object An object containing an edgelist
#' @param components Components to compute A node projection for
#' @param verbose Print messages
#' @param ... Not yet implemented
#'
#' @rdname graph-conversion
#'
#' @return An A-node-projected graph
#'
#' @export
#'
edgelist_to_simple_Anode_graph <- function (
  object,
  ...
) {
  UseMethod(generic = "edgelist_to_simple_Anode_graph", object = object)
}


#' Differential analysis (polarization)
#'
#' Runs differential analysis on Moran's Z polarization scores generated with the
#' \code{Pixelator} data processing pipeline.
#'
#' If you are working with a \code{Seurat} object containing a \code{\link{CellGraphAssay}},
#' the polarization scores are accessed directly from the \code{\link{CellGraphAssay}}.
#' A character vector or factor must be selected with \code{contrast_column} from the
#' input data (or @meta.data slot from a \code{Seurat} object) which holds the groups
#' to run the test for. The \code{target} and \code{reference} parameters should refer
#' to the names of the two groups used for the comparison and these names should be present in
#' the \code{contrast_column}.
#'
#' @section Additional groups:
#' The test is always computed between \code{target} and \code{reference}, but it is possible
#' to add additional grouping variables with \code{group_vars}. If \code{group_vars} is used,
#' the test will be computed within each combination of groups. For instance, if we have annotated
#' cells into cell type populations across two conditions defined by \code{target} and \code{reference},
#' we can pass the name of a cell annotation column with \code{group_vars} to run the test
#' for each cell type.
#'
#' @concept DA
#' @family DA-methods
#'
#' @param object An object containing polarization scores
#' @param target The name of the target group
#' @param reference The name of the reference group
#' @param contrast_column The name of the column where the group labels are stored.
#' This column must include \code{target} and \code{reference}.
#' @param group_vars An optional character vector with column names to split the tests by.
#' @param alternative One of 'two.sided', 'less' or 'greater' (see \code{?wilcox.test} for details)
#' @param conf_int Should confidence intervals be computed? (see \code{?wilcox.test} for details)
#' @param p_adjust_method One of "bonferroni", "holm", "hochberg", "hommel", "BH", "BY" or "fdr".
#' (see \code{?p.adjust} for details)
#' @param verbose Print messages
#' @param ... Not yet implemented
#'
#' @rdname RunDPA
#'
#' @return A \code{tbl_df} object with test results
#'
#' @export
#'
RunDPA <- function (
  object,
  ...
) {
  UseMethod(generic = "RunDPA", object = object)
}


#' Differential analysis (colocalization)
#'
#' Runs differential analysis on Moran's Z colocalization scores generated with the
#' \code{Pixelator} data processing pipeline.
#'
#' If you are working with a \code{Seurat} object containing a \code{\link{CellGraphAssay}},
#' the polarization scores are accessed directly from the \code{\link{CellGraphAssay}}.
#' A character vector or factor must be selected with \code{contrast_column} from the
#' input data (or @meta.data slot from a \code{Seurat} object) which holds the groups
#' to run the test for. The \code{target} and \code{reference} parameters should refer
#' to the names of the two groups used for the comparison and these names should be present in
#' the \code{contrast_column}.
#'
#' @section Additional groups:
#' The test is always computed between \code{target} and \code{reference}, but it is possible
#' to add additional grouping variables with \code{group_vars}. If \code{group_vars} is used,
#' the test will be computed within each combination of groups. For instance, if we have annotated
#' cells into cell type populations across two conditions defined by \code{target} and \code{reference},
#' we can pass the name of a cell annotation column with \code{group_vars} to run the test
#' for each cell type.
#'
#' @concept DA
#' @family DA-methods
#'
#' @param object An object containing colocalization scores
#' @param target The name of the target group
#' @param reference The name of the reference group
#' @param contrast_column The name of the column where the group labels are stored.
#' This column must include \code{target} and \code{reference}.
#' @param group_vars An optional character vector with column names to group the tests by.
#' @param alternative One of 'two.sided', 'less' or 'greater' (see \code{?wilcox.test} for details)
#' @param conf_int Should confidence intervals be computed? (see \code{?wilcox.test} for details)
#' @param p_adjust_method One of "bonferroni", "holm", "hochberg", "hommel", "BH", "BY" or "fdr".
#' (see \code{?p.adjust} for details)
#' @param verbose Print messages
#' @param ... Not yet implemented
#'
#' @rdname RunDCA
#'
#' @return A \code{tbl_df} object with test results
#'
#' @export
#'
RunDCA <- function (
  object,
  ...
) {
  UseMethod(generic = "RunDCA", object = object)
}


#' Convert polarization score table to an Assay or Assay5
#'
#' @section Behavior:
#'
#' Takes an object with polarization scores in long format and returns an
#' object with polarization scores in a wide format. The polarization score
#' table includes Moran's I and Z scores along with p-values for each marker
#' and component.
#'
#' The wide format is an array-like object with dimensions markers x components,
#' where each cell is filled with a polarization score. Scores that are missing
#' from the polarization score table are replaced with 0's.
#'
#' Different outputs are returned depending on the input object type:
#'
#' \itemize{
#'    \item{
#'      \code{tibble/data.frame}: returns a matrix with markers in rows and components
#'      in columns
#'    }
#'    \item{
#'      \code{CelGraphAssay}: returns an Assay with markers in rows and
#'      components in columns
#'    }
#'    \item{
#'      \code{Seurat} object: returns the Seurat object with a new Assay
#'      with markers in rows and components in columns
#'    }
#' }
#'
#' As many functions provided in Seurat works on \code{Assay} objects, it is
#' sometimes convenient to make this conversion. For instance, if we want to
#' compute a UMAP on the polarization scores with \code{RunUMAP}, we need the
#' values to be formatted in an \code{Assay}. This also makes it possible
#' to use various visualization functions such as \code{VlnPlot} or \code{FeaturePlor}
#' to show the distribution of polarization scores.
#'
#' @param object An object with polarization scores
#' @param values_from What column to pick polarization scores from. Either
#' "morans_i" or "morans_z"
#' @param ... Not yet implemented
#'
#' @rdname PolarizationScoresToAssay
#' @family Spatial metrics conversion methods
#'
#' @export
#'
PolarizationScoresToAssay <- function (
  object,
  ...
) {
  UseMethod(generic = "PolarizationScoresToAssay", object = object)
}


#' Convert colocalization score table to an Assay or Assay5
#'
#' @section Behavior:
#'
#' Takes an object with colocalization scores in long format and returns an
#' object with colocalization scores in a wide format. The colocalization score
#' table includes various colocalization scores along with p-values for each pair
#' markers and component.
#'
#' The wide format is an array-like object with dimensions (markers_1 * marker_2) x components,
#' where each cell is filled with a polarization score. Scores that are missing
#' from the colocalization score table are replaced with 0's.
#'
#' Different outputs are returned depending on the input object type:
#'
#' \itemize{
#'    \item{
#'      \code{tibble/data.frame}: returns a matrix with marker pairs in rows and components
#'      in columns
#'    }
#'    \item{
#'      \code{CelGraphAssay}: returns an Assay with marker pairs in rows and
#'      components in columns
#'    }
#'    \item{
#'      \code{Seurat} object: returns the Seurat object with a new Assay
#'      with marker pairs in rows and components in columns
#'    }
#' }
#'
#' As many functions provided in Seurat works on \code{Assay} objects, it is
#' sometimes convenient to make this conversion. For instance, if we want to
#' compute a UMAP on the colocalization scores with \code{RunUMAP}, we need the
#' values to be formatted in an \code{Assay}. This also makes it possible
#' to use various visualization functions such as \code{VlnPlot} or \code{FeaturePlor}
#' to show the distribution of colocalization scores.
#'
#' @param object An object with colocalization scores
#' @param values_from What column to pick colocalization scores from. One of
#' "pearson" or "pearson_z"
#' @param ... Not yet implemented
#'
#' @rdname ColocalizationScoresToAssay
#' @family Spatial metrics conversion methods
#'
#' @export
#'
ColocalizationScoresToAssay <- function (
  object,
  ...
) {
  UseMethod(generic = "ColocalizationScoresToAssay", object = object)
}

#' PolarizationScores
#'
#' Get and set polarization scores for a \code{\link{CellGraphAssay}},
#' \code{\link{CellGraphAssay5}} or a \code{Seurat} object
#'
#' @param object An object with polarization scores
#' @param ... Not implemented
#'
#' @rdname PolarizationScores
#' @family spatial metrics
#'
#' @return \code{PolarizationScores}: Polarization scores
#'
#' @export
#'
PolarizationScores <- function (
  object,
  ...
) {
  UseMethod(generic = "PolarizationScores", object = object)
}


#' @param value A \code{tbl_df} with polarization scores
#'
#' @rdname PolarizationScores
#'
#' @return \code{PolarizationScores<-}: An object with polarization scores
#' updated
#'
#' @export
#'
"PolarizationScores<-" <- function (
  object,
  ...,
  value
) {
  UseMethod(generic = "PolarizationScores<-", object = object)
}


#' ColocalizationScores
#'
#' Get and set colocalization scores for a \code{\link{CellGraphAssay}},
#' \code{\link{CellGraphAssay5}} or a \code{Seurat} object
#'
#' @param object An object with polarization scores
#' @param ... Not implemented
#'
#' @rdname ColocalizationScores
#' @family spatial metrics
#'
#' @return \code{ColocalizationScores}: Colocalization scores
#'
#' @export
#'
ColocalizationScores <- function (
  object,
  ...
) {
  UseMethod(generic = "ColocalizationScores", object = object)
}


#' @param value A \code{tbl_df} with colocalization scores
#'
#' @rdname ColocalizationScores
#'
#' @return \code{ColocalizationScores<-}: An object with colocalization scores
#' updated
#'
#' @export
#'
"ColocalizationScores<-" <- function (
  object,
  ...,
  value
) {
  UseMethod(generic = "ColocalizationScores<-", object = object)
}


#' FS map
#'
#' Get and set \code{fs_map} tibble for a \code{\link{CellGraphAssay}},
#' \code{\link{CellGraphAssay5}} or a \code{Seurat} object
#'
#' @param object An object containing arrow data
#' @param ... Not implemented
#'
#' @rdname FSMap
#'
#' @seealso [CellGraphs()] for getting/setting \code{\link{CellGraph}} lists
#' and [PolarizationScores()],[ColocalizationScores()] for getting/setting
#' spatial metrics
#'
#' @return \code{FSMap}: An \code{fs_map} tibble
#'
#' @export
#'
FSMap <- function (
  object,
  ...
) {
  UseMethod(generic = "FSMap", object = object)
}


#' @param value A \code{tbl_df}
#'
#' @rdname FSMap
#'
#' @return \code{FSMap<-}: An object with an updated \code{fs_map} tibble
#'
#' @export
#'
"FSMap<-" <- function (
  object,
  ...,
  value
) {
  UseMethod(generic = "FSMap<-", object = object)
}
