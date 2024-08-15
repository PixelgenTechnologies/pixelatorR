#' @include generics.R
NULL

#' Convert objects to a \code{\link{CellGraphAssay5}}
#'
#' @param x An object to convert to class \code{\link{CellGraphAssay5}}
#' @param ... Arguments passed to other methods
#' @rdname as.CellGraphAssay5
#' @export as.CellGraphAssay5
as.CellGraphAssay5 <- function(
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
as.CellGraphAssay <- function(
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
CellGraphs <- function(
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
"CellGraphs<-" <- function(
  object,
  ...,
  value
) {
  UseMethod(generic = "CellGraphs<-", object = object)
}

#' Load CellGraphs
#'
#' Loads a list of \code{\link{CellGraph}} objects.
#'
#' Graphs can be loaded as one of 'bipartite', 'Anode' or 'linegraph'. See details about each of
#' these graph representations in the sections below.
#'
#' \code{LoadCellGraphs} will only work if the PXL file path(s) stored in the object are valid.
#' You can check the PXL file paths stored in a \code{Seurat}, a \code{CellGraphAssay} or a
#' \code{CellGraphAssay5} object with the \code{\link{FSMap}} method. If the PXL file(s) are
#' invalid, an error will be thrown. Visit \code{\link{RestorePaths}} for more information on how
#' to update the PXL file paths.
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
LoadCellGraphs <- function(
  object,
  ...
) {
  UseMethod(generic = "LoadCellGraphs", object = object)
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
RemoveCellGraphs <- function(
  object,
  ...
) {
  UseMethod(generic = "RemoveCellGraphs", object = object)
}

#' Edge Rank Plot
#'
#' Plots the number of edges/molecules per component against the component rank
#'
#' @param object A \code{data.frame}-like object or a \code{Seurat} object
#' @param ... Parameters passed to other methods
#'
#' @rdname MoleculeRankPlot
#'
#' @return A \code{ggplot} object
#'
#' @export MoleculeRankPlot
#'
MoleculeRankPlot <- function(
  object,
  ...
) {
  UseMethod(generic = "MoleculeRankPlot", object = object)
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
#' @param as_frequency Plot frequencies instead of counts
#' @param stack Create a stacked bar plot. Only has an effect if a \code{group_by}
#' variable is provided.
#' @param ... Not yet implemented
#'
#' @rdname CellCountPlot
#'
#' @return A \code{ggplot} object
#'
#' @export CellCountPlot
#'
CellCountPlot <- function(
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
TauPlot <- function(
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
ComputeLayout <- function(
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
KeepLargestComponent <- function(
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
edgelist_to_simple_Anode_graph <- function(
  object,
  ...
) {
  UseMethod(generic = "edgelist_to_simple_Anode_graph", object = object)
}


#' Differential analysis (polarity)
#'
#' Runs differential analysis on polarity scores generated with the
#' \code{Pixelator} data processing pipeline.
#'
#' If you are working with a \code{Seurat} object created with pixelatorR that contains
#' a \code{\link{CellGraphAssay}}, the polarity scores are accessed directly from
#' the \code{\link{CellGraphAssay}} (see \code{\link{PolarizationScores}}).
#'
#' The input object should contain a \code{contrast_column} (character vector or factor)
#' that includes information about the groups to compare. A typical example is a column
#' with sample labels, for instance: "control", "stimulated1", "stimulated2". If the input
#' object is a \code{Seurat} object, the \code{contrast_column} should be available in
#' the \code{meta.data} slot. For those familiar with \code{FindMarkers} from Seurat,
#' \code{contrast_column} is equivalent to the \code{group.by} parameter.
#'
#' The \code{targets} parameter specifies a character vector with the names of the groups
#' to compare \code{reference}. \code{targets} can be a single group name or a vector of
#' group names while \code{reference} can only refer to a single group. Both \code{targets}
#' and \code{reference} should be present in the \code{contrast_column}. These parameters
#' are similar to the \code{ident.1} and \code{ident.2} parameters in \code{FindMarkers}.
#'
#' @section Additional groups:
#' The test is always computed between \code{targets} and \code{reference}, but it is possible
#' to add additional grouping variables with \code{group_vars}. If \code{group_vars} is used,
#' each comparison is split into groups defined by the \code{group_vars}. For instance, if we
#' have annotated cells into cell type populations and saved these annotations in a \code{meta.data}
#' column called "cell_type", we can pass "cell_type" to \code{group_vars="cell_type"} to split
#' tests across each cell type.
#'
#' @section Types of comparisons:
#' Consider a scenario where we have a Seurat object (\code{seurat_object}) with MPX data.
#' \code{seurat_object} contains a \code{meta.data} column called "sampleID" that holds
#' information about what samples the MPX components originated from. This column could have
#' three sample IDs: "control", "stimulated1" and "stimulated2". In addition, we have a column
#' called "cell_type" that holds information about the cell type identity of each MPX component.
#'
#' 1. If we want to compare the "stimulated1" group to the "control" group:
#' \preformatted{
#' dpa_markers <- RunDPA(object = seurat_object,
#'                       contrast_column = "sampleID",
#'                       reference = "control",
#'                       targets = "stimulated1")
#' }
#'
#' 2. If we want to compare the "stimulated1" and "stimulated2" groups to the "control" group:
#' \preformatted{
#' dpa_markers <- RunDPA(object = seurat_object,
#'                       contrast_column = "sampleID",
#'                       reference = "control",
#'                       targets = c("stimulated1", "stimulated2"))
#' }
#'
#' 3. If we want to compare the "stimulated1" and "stimulated2" groups to the "control" group, and split
#' the tests by cell type:
#' \preformatted{
#' dpa_markers <- RunDPA(object = seurat_object,
#'                      contrast_column = "sampleID",
#'                      reference = "control",
#'                      targets = c("stimulated1", "stimulated2"),
#'                      group_vars = "cell_type")
#' }
#'
#' @concept DA
#' @family DA-methods
#'
#' @param object An object containing polarity scores
#' @param contrast_column The name of the column where the group labels are stored.
#' This column must include \code{target} and \code{reference}.
#' @param targets The name of the target groups. These groups will be compared to the reference group.
#' If the value is set to \code{NULL} (default), all groups available in \code{contrast_column} except
#' \code{reference} will be compared to the \code{reference} group.
#' @param reference The name of the reference group
#' @param group_vars An optional character vector with column names to group the tests by.
#' @param polarity_metric The polarity metric to use. Currently, you can select one of "morans_z" (default)
#' or "morans_i".
#' @param min_n_obs Minimum number of observations allowed in a group. Target groups with less
#' observations than \code{min_n_obs} will be skipped.
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
RunDPA <- function(
  object,
  ...
) {
  UseMethod(generic = "RunDPA", object = object)
}


#' Differential analysis (colocalization)
#'
#' Runs differential analysis on colocalization scores generated with the
#' \code{Pixelator} data processing pipeline.
#'
#' If you are working with a \code{Seurat} object created with pixelatorR that contains
#' a \code{\link{CellGraphAssay}}, the colocalization scores are accessed directly from
#' the \code{\link{CellGraphAssay}} (see \code{\link{ColocalizationScores}}).
#'
#' The input object should contain a \code{contrast_column} (character vector or factor)
#' that includes information about the groups to compare. A typical example is a column
#' with sample labels, for instance: "control", "stimulated1", "stimulated2". If the input
#' object is a \code{Seurat} object, the \code{contrast_column} should be available in
#' the \code{meta.data} slot. For those familiar with \code{FindMarkers} from Seurat,
#' \code{contrast_column} is equivalent to the \code{group.by} parameter.
#'
#' The \code{targets} parameter specifies a character vector with the names of the groups
#' to compare \code{reference}. \code{targets} can be a single group name or a vector of
#' group names while \code{reference} can only refer to a single group. Both \code{targets}
#' and \code{reference} should be present in the \code{contrast_column}. These parameters
#' are similar to the \code{ident.1} and \code{ident.2} parameters in \code{FindMarkers}.
#'
#' @section Additional groups:
#' The test is always computed between \code{targets} and \code{reference}, but it is possible
#' to add additional grouping variables with \code{group_vars}. If \code{group_vars} is used,
#' each comparison is split into groups defined by the \code{group_vars}. For instance, if we
#' have annotated cells into cell type populations and saved these annotations in a \code{meta.data}
#' column called "cell_type", we can pass "cell_type" to \code{group_vars="cell_type"} to split
#' tests across each cell type.
#'
#' @section Types of comparisons:
#' Consider a scenario where we have a Seurat object (\code{seurat_object}) with MPX data.
#' \code{seurat_object} contains a \code{meta.data} column called "sampleID" that holds
#' information about what samples the MPX components originated from. This column could have
#' three sample IDs: "control", "stimulated1" and "stimulated2". In addition, we have a column
#' called "cell_type" that holds information about the cell type identity of each MPX component.
#'
#' 1. If we want to compare the "stimulated1" group to the "control" group:
#' \preformatted{
#' dca_markers <- RunDCA(object = seurat_object,
#'                       contrast_column = "sampleID",
#'                       reference = "control",
#'                       targets = "stimulated1")
#' }
#'
#' 2. If we want to compare the "stimulated1" and "stimulated2" groups to the "control" group:
#' \preformatted{
#' dca_markers <- RunDCA(object = seurat_object,
#'                       contrast_column = "sampleID",
#'                       reference = "control",
#'                       targets = c("stimulated1", "stimulated2"))
#' }
#'
#' 3. If we want to compare the "stimulated1" and "stimulated2" groups to the "control" group, and split
#' the tests by cell type:
#' \preformatted{
#' dca_markers <- RunDCA(object = seurat_object,
#'                      contrast_column = "sampleID",
#'                      reference = "control",
#'                      targets = c("stimulated1", "stimulated2"),
#'                      group_vars = "cell_type")
#' }
#'
#' @concept DA
#' @family DA-methods
#'
#' @param object An object containing colocalization scores
#' @param contrast_column The name of the column where the group labels are stored.
#' This column must include \code{target} and \code{reference}.
#' @param targets The name of the target groups. These groups will be compared to the reference group.
#' If the value is set to \code{NULL} (default), all groups available in \code{contrast_column} except
#' \code{reference} will be compared to the \code{reference} group.
#' @param reference The name of the reference group
#' @param group_vars An optional character vector with column names to group the tests by.
#' @param coloc_metric The colocalization metric to use. Currently, you can select one of "pearson_z" (default)
#' or "pearson".
#' @param min_n_obs Minimum number of observations allowed in a group. Target groups with less
#' observations than \code{min_n_obs} will be skipped.
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
RunDCA <- function(
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
PolarizationScoresToAssay <- function(
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
ColocalizationScoresToAssay <- function(
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
PolarizationScores <- function(
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
"PolarizationScores<-" <- function(
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
ColocalizationScores <- function(
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
"ColocalizationScores<-" <- function(
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
FSMap <- function(
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
"FSMap<-" <- function(
  object,
  ...,
  value) {
  UseMethod(generic = "FSMap<-", object = object)
}


#' Normalize MPX data
#'
#' Normalizes MPX data using the specified method. The normalization method can be one of "dsb" or "CLR".
#'
#' CLR can be used to normalize MPX data using the centered log-ratio transformation in which
#' the assumption is that the geometric mean of the marker abundance is constant across cells
#' (e.g cell lines). This assumption might not hold for datasets from samples from different
#' sources or having a variable cell type composition. In addition, CLR does not take the noise
#' from unspecific binding of antibodies into account.
#'
#' For these reasons, the dsb normalization method can be a useful alternative in mixed-population
#' datasets. dsb normalizes marker counts based on their abundance in a negative population across
#' all cells and regresses out a per-cell noise component based on isotype controls and non-specific
#' marker abundance.
#'
#' @references Mulè, M.P., Martins, A.J. & Tsang, J.S. Normalizing and denoising protein expression
#' data from droplet-based single cell profiling. Nat Commun 13, 2099 (2022).
#' \url{https://doi.org/10.1038/s41467-022-29356-8}
#'
#' @param object An object.
#' @param method The normalization method to use. Can be "dsb" or "clr".
#' @param isotype_controls A vector of isotype controls to use for normalization.
#' @param assay Name of assay to use; defaults to the default assay.
#' @param ... Additional arguments. Currently not used.
#'
#' @rdname NormalizeMPX
#'
#' @return An object with normalized MPX data.
#'
#' @export
#'
NormalizeMPX <- function(
  object,
  method = c("dsb", "clr"),
  isotype_controls = c("mIgG1", "mIgG2a", "mIgG2b"),
  assay = NULL,
  ...
) {
  UseMethod(generic = "NormalizeMPX", object = object)
}


#' Restore PXL file paths
#'
#' Updates the PXL file paths in an object with the paths in the specified directory.
#'
#' @details
#' Seurat objects created with \code{pixelatorR} store the MPX data
#' in a \code{\link{CellGraphAssay}} or \code{\link{CellGraphAssay5}} object.
#' For some analytical tasks, such as component visualization, you need
#' to load the component graphs in memory with \code{\link{LoadCellGraphs}}.
#'
#' \code{\link{LoadCellGraphs}} reads data from the edgelist(s) stored in the PXL
#' files(s) that the Seurat object was created from (see \code{\link{ReadMPX_Seurat}}).
#' This means that the PXL files must be accessible in the \code{\link{CellGraphAssay}}
#' or \code{\link{CellGraphAssay5}} when you load the component graphs.
#' You can check the PXL file paths with \code{\link{FSMap}}.
#'
#' If the PXL files are moved or copied to a different directory, \code{\link{LoadCellGraphs}}
#' will fail because the PXL file paths are invalid. \code{RestorePaths} can
#' be useful in these situations to update the PXL file paths in the Seurat object.
#' All you need is to provide the path to the directory where the PXL files are located.
#'
#' @section Exported Seurat objects:
#' A typical situation when this function is useful is when you export a Seurat object
#' created with \code{pixelatorR} to an RDS file, and share it with someone else. To ensure
#' that all data is available for the recipient, you also need to share the raw PXL files
#' and make sure that the recipient updates the PXL file paths in the Seurat object with
#' \code{RestorePaths}.
#'
#' For instance, let's assume that we have a Seurat object \code{seurat_obj} created with
#' \code{pixelatorR} and that the PXL files are stored in the directory \code{DATA_DIR}
#' on your local system. You can export the Seurat object to an RDS file with:
#' \preformatted{
#' saveRDS(seurat_obj, file = file.path(DATA_DIR, "seurat_obj.rds"))
#' }
#' The content of DATA_DIR might look like this:
#' \preformatted{
#' DATA_DIR
#' ├── seurat_obj.rds
#' └── Sample01_pbmc.analysis.pxl
#' }
#' Now you can share the entire \code{DATA_DIR} folder with someone else. When the recipient
#' loads the RDS file in R, the PXL file paths in the Seurat object need to be updated,
#' because recipients \code{DATA_DIR} will be different from yours. Assuming that the recipient
#' has stored the PXL file and the RDS object in the directory \code{RECIPIENT_DATA_DIR},
#' the recipient can update the PXL file paths with:
#' \preformatted{
#' seurat_obj <- readRDS(file = file.path(RECIPIENT_DATA_DIR, "seurat_obj.rds"))
#' seurat_obj <- RestorePaths(seurat_obj, pxl_files_dir = RECIPIENT_DATA_DIR)
#' }
#'
#' @param object An object
#'
#' @rdname RestorePaths
#'
#' @return An object with updated PXL file paths
#'
#' @export
#'
RestorePaths <- function(
  object,
  ...
) {
  UseMethod(generic = "RestorePaths", object = object)
}
