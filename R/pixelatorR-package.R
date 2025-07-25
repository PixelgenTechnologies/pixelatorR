#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom arrow arrow_table
#' @importFrom arrow open_dataset
#' @importFrom arrow read_parquet
#' @importFrom arrow schema
#' @importFrom arrow to_duckdb
#' @importFrom arrow write_dataset
#' @importFrom cli cat_line
#' @importFrom cli cli_alert
#' @importFrom cli cli_alert_danger
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_alert_success
#' @importFrom cli cli_alert_warning
#' @importFrom cli cli_div
#' @importFrom cli cli_end
#' @importFrom cli cli_h2
#' @importFrom cli cli_progress_bar
#' @importFrom cli cli_progress_done
#' @importFrom cli cli_progress_update
#' @importFrom cli cli_rule
#' @importFrom cli cli_status
#' @importFrom cli cli_status_clear
#' @importFrom cli cli_status_update
#' @importFrom cli cli_ul
#' @importFrom cli cli_warn
#' @importFrom cli col_blue
#' @importFrom cli col_br_blue
#' @importFrom cli col_br_magenta
#' @importFrom cli col_br_red
#' @importFrom cli col_green
#' @importFrom cli col_red
#' @importFrom cli symbol
#' @importFrom dplyr %>%
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr any_of
#' @importFrom dplyr arrange
#' @importFrom dplyr between
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom dplyr case_when
#' @importFrom dplyr collect
#' @importFrom dplyr compute
#' @importFrom dplyr contains
#' @importFrom dplyr copy_to
#' @importFrom dplyr count
#' @importFrom dplyr cross_join
#' @importFrom dplyr cur_group_id
#' @importFrom dplyr desc
#' @importFrom dplyr distinct
#' @importFrom dplyr do
#' @importFrom dplyr everything
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr group_by_all
#' @importFrom dplyr group_by_at
#' @importFrom dplyr group_data
#' @importFrom dplyr group_keys
#' @importFrom dplyr group_split
#' @importFrom dplyr if_else
#' @importFrom dplyr inner_join
#' @importFrom dplyr is.grouped_df
#' @importFrom dplyr left_join
#' @importFrom dplyr matches
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_all
#' @importFrom dplyr mutate_if
#' @importFrom dplyr n
#' @importFrom dplyr n_distinct
#' @importFrom dplyr pick
#' @importFrom dplyr pull
#' @importFrom dplyr reframe
#' @importFrom dplyr relocate
#' @importFrom dplyr rename
#' @importFrom dplyr rename_with
#' @importFrom dplyr right_join
#' @importFrom dplyr row_number
#' @importFrom dplyr rowwise
#' @importFrom dplyr select
#' @importFrom dplyr slice
#' @importFrom dplyr slice_head
#' @importFrom dplyr summarise
#' @importFrom dplyr summarize
#' @importFrom dplyr tbl
#' @importFrom dplyr ungroup
#' @importFrom dplyr vars
#' @importFrom dplyr where
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 expansion
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 layer_scales
#' @importFrom ggplot2 position_dodge
#' @importFrom ggplot2 position_stack
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_color_identity
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_color_viridis_c
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_identity
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_size
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_y_log10
#' @importFrom ggplot2 sym
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme_void
#' @importFrom ggraph geom_edge_link
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph ggraph
#' @importFrom glue glue
#' @importFrom glue glue_sql
#' @importFrom graphics par
#' @importFrom grDevices as.raster
#' @importFrom grDevices col2rgb
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices dev.off
#' @importFrom grid gpar
#' @importFrom grid textGrob
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph connect
#' @importFrom igraph gsize
#' @importFrom igraph vertex_attr_names
#' @importFrom Matrix colSums
#' @importFrom Matrix rowSums
#' @importFrom methods as
#' @importFrom methods getMethod
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom methods setClass
#' @importFrom methods setClassUnion
#' @importFrom methods setMethod
#' @importFrom methods show
#' @importFrom methods slot
#' @importFrom methods slot<-
#' @importFrom methods slotNames
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom patchwork plot_annotation
#' @importFrom patchwork plot_layout
#' @importFrom patchwork plot_spacer
#' @importFrom patchwork wrap_elements
#' @importFrom patchwork wrap_plots
#' @importFrom pbapply pblapply
#' @importFrom R6 R6Class
#' @importFrom rlang on_load
#' @importFrom rlang run_on_load
#' @importFrom SeuratObject Cells
#' @importFrom SeuratObject CreateAssay5Object
#' @importFrom SeuratObject CreateAssayObject
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom SeuratObject DefaultAssay
#' @importFrom SeuratObject FetchData
#' @importFrom SeuratObject JoinLayers
#' @importFrom SeuratObject Key
#' @importFrom SeuratObject Key<-
#' @importFrom SeuratObject LayerData
#' @importFrom SeuratObject LayerData<-
#' @importFrom SeuratObject Layers
#' @importFrom SeuratObject RenameCells
#' @importFrom SeuratObject VariableFeatures
#' @importFrom SeuratObject VariableFeatures<-
#' @importFrom stats as.formula
#' @importFrom stats cor
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats p.adjust
#' @importFrom stats pnorm
#' @importFrom stats prcomp
#' @importFrom stats quantile
#' @importFrom stats reformulate
#' @importFrom stats rnorm
#' @importFrom stats wilcox.test
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_sub
#' @importFrom tibble as_tibble
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom tibble tibble
#' @importFrom tidygraph %E>%
#' @importFrom tidygraph %N>%
#' @importFrom tidygraph as_tbl_graph
#' @importFrom tidygraph bind_edges
#' @importFrom tidygraph to_components
#' @importFrom tidyr crossing
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
#' @importFrom tidyr unite
#' @importFrom utils capture.output
#' @importFrom utils head
#' @importFrom utils menu
#' @importFrom utils unzip
## usethis namespace: end
NULL
