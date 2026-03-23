# Package index

## Function reference

### Load data

- [`ReadMPX_Seurat()`](ReadMPX_Seurat.md) :

  Load data from PXL file into a `Seurat` object

- [`ReadMPX_counts()`](ReadMPX_counts.md) : Read a count matrix from a
  pxl file

- [`ReadMPX_item()`](ReadMPX_item.md)
  [`ReadMPX_polarization()`](ReadMPX_item.md)
  [`ReadMPX_colocalization()`](ReadMPX_item.md)
  [`ReadMPX_edgelist()`](ReadMPX_item.md) : Read an MPX data item

- [`ReadMPX_arrow_edgelist()`](ReadMPX_arrow_edgelist.md) : Read
  edgelists from a PXL file containing MPX data

- [`ReadMPX_layouts()`](ReadMPX_layouts.md) : Load layouts from an PXL
  file containing MPX data

- [`ReadMPX_metadata()`](ReadMPX_metadata.md) : Read metadata from a PXL
  file

- [`ReadPNA_Seurat()`](ReadPNA_Seurat.md) :

  Load data from PNA PXL file into a `Seurat` object

- [`ReadPNA_counts()`](ReadPNA_counts.md) : Read a count matrix from a
  PXL file with PNA data

- [`ReadPNA_edgelist()`](ReadPNA_edgelist.md) : Load the edge list from
  a PNA PXL file

- [`ReadPNA_proximity()`](ReadPNA_proximity.md) : Load the Proximity
  scores table from a PNA PXL file

- [`ReadPNA_layouts()`](ReadPNA_layouts.md) : Load layouts from a PNA
  PXL file

- [`ReadPNA_metadata()`](ReadPNA_metadata.md) : Read metadata from a PNA
  PXL file

### PXL file database

- [`PixelDB`](PixelDB.md) : PXL database class

### Visualization

- [`CellCountPlot()`](CellCountPlot.md) : Plot cell counts per group
- [`MoleculeRankPlot()`](MoleculeRankPlot.md) : Edge Rank Plot
- [`EdgeRankPlot()`](EdgeRankPlot.md) **\[deprecated\]** : Edge Rank
  Plot
- [`TauPlot()`](TauPlot.md) : Plot UMIs per UPIa for quality control
- [`DensityScatterPlot()`](DensityScatterPlot.md) : Create a density
  scatter / pseudocolor plot.
- [`AbundanceColocalizationPlot()`](AbundanceColocalizationPlot.md) :
  Create an abundance/colocalization scatterplot
- [`Plot2DGraph()`](Plot2DGraph.md) : Plot 2D graph layouts
- [`Plot2DGraphM()`](Plot2DGraphM.md) : Plot multiple markers on
  multiple graphs
- [`Plot3DGraph()`](Plot3DGraph.md) : Plot 3D graph layouts

### Graph layouts

- [`ComputeLayout()`](ComputeLayout.md) : Compute a graph layout
- [`layout_with_weighted_pmds()`](layout_with_weighted_pmds.md) :
  Weighted PMDS
- [`fast_pmds()`](fast_pmds.md) : Fast pMDS implementation using
  RSpectra
- [`center_layout_coordinates()`](layout_coordinates_utils.md)
  [`normalize_layout_coordinates()`](layout_coordinates_utils.md)
  [`project_layout_coordinates_on_unit_sphere()`](layout_coordinates_utils.md)
  : Layout Coordinates Utility Functions
- [`render_rotating_layout()`](render_rotating_layout.md)
  **\[experimental\]** : Create a rotating 3D layout video
- [`scale_layout()`](scale_layout.md) : Scale layout coordinates

### Data processing

- [`NormalizeMPX()`](NormalizeMPX.md) **\[superseded\]** : Normalize MPX
  data
- [`KeepLargestComponent()`](KeepLargestComponent.md) : Keep largest
  component
- [`RemoveCellGraphs()`](RemoveCellGraphs.md) : Remove CellGraphs
- [`LoadCellGraphs()`](LoadCellGraphs.md) : Load CellGraphs
- [`cos_dist()`](cos_dist.md) : Calculate cosine distances between two
  sets of coordinates
- [`pack_2bits()`](two_bit_encoding.md)
  [`unpack_2bits()`](two_bit_encoding.md) : Encode/decode DNA sequences

### Differential analysis

- [`RunDAA()`](RunDAA.md) : Differential analysis (abundance)
- [`RunDPA()`](RunDPA.md) : Differential analysis (polarity)
- [`RunDCA()`](RunDCA.md) : Differential analysis (colocalization)
- [`DifferentialProximityAnalysis()`](DifferentialProximityAnalysis.md)
  : Differential analysis of proximity scores
- [`ColocalizationHeatmap()`](ColocalizationHeatmap.md) : Plot a
  colocalization heatmap

### Proximity scores

- [`PolarizationScoresToAssay()`](PolarizationScoresToAssay.md) :
  Convert polarization score table to an Assay or Assay5
- [`ColocalizationScoresToAssay()`](ColocalizationScoresToAssay.md) :
  Convert colocalization score table to an Assay or Assay5
- [`ProximityScoresToAssay()`](ProximityScoresToAssay.md) : Convert
  proximity score table to an Assay or Assay5
- [`FilterProximityScores()`](FilterProximityScores.md) : Filter
  proximity scores
- [`SummarizeProximityScores()`](SummarizeProximityScores.md) :
  Summarize proximity scores

### Sequencing saturation

- [`sequencing_saturation()`](sequencing_saturation.md) : Calculate
  Sequencing Saturation
- [`SequenceSaturationCurve()`](SequenceSaturationCurve.md) : Simulate
  Sequencing Saturation Curve
- [`approximate_edge_saturation()`](approximate_edge_saturation.md) :
  Compute approximate edge saturation
- [`approximate_node_saturation()`](approximate_node_saturation.md) :
  Compute approximate node saturation
- [`approximate_saturation_curve()`](approximate_saturation_curve.md) :
  Compute approximate saturation curve
- [`downsample_to_parquet()`](downsample_to_parquet.md) : Filter an
  edgelist by downsampling read counts
- [`lcc_sizes()`](lcc_sizes.md) : Calculate the sizes of the LCCs from
  an edgelist
- [`lcc_curve()`](lcc_curve.md) **\[experimental\]** : Compute LCC sizes
  for downsampled edgelists

### Cell Annotation

- [`AnnotateCells()`](AnnotateCells.md) : Automatic annotation of cell
  types
- [`read_pbmc_reference()`](read_pbmc_reference.md) : Read PBMC cell
  annotation reference

### Color themes

- [`color_discrete_pixelgen()`](color-themes.md)
  [`color_sequential_pixelgen()`](color-themes.md)
  [`color_divergent_pixelgen()`](color-themes.md) : Pixelgen ggplot2
  color themes
- [`PixelgenAccentColors()`](PixelgenAccentColors.md) : Get Pixelgen
  accent colors
- [`PixelgenGradient()`](PixelgenGradient.md) : Get Pixelgen gradient
  colors
- [`PixelgenPalette()`](PixelgenPalette.md) : Get Pixelgen palette
  colors
- [`Pixelgen_accent_colors`](Pixelgen_accent_colors.md) : Pixelgen
  accent colors
- [`show_accent_colors()`](show_accent_colors.md) : Show accent colors
- [`theme_pixelgen()`](theme_pixelgen.md) : Pixelgen theme
- [`create_discrete_palette()`](create_discrete_palette.md) : Create
  sample palette
- [`Pixelgen_cell_palette`](Pixelgen_cell_palette.md) : Pixelgen cell
  type colors

### Doublet detection

- [`FindAnnoyNeighbors()`](FindAnnoyNeighbors.md) : Approximate nearest
  neighbors using Annoy
- [`SimulateDoublets()`](SimulateDoublets.md) : Simulate doublets
- [`PredictDoublets()`](PredictDoublets.md) **\[experimental\]** :
  Predict doublets in a Seurat object

### Patch analysis

- [`identify_markers_for_patch_analysis()`](identify_markers_for_patch_analysis.md)
  : Identify population-specific markers for patch detection
- [`patch_detection()`](patch_detection.md) **\[experimental\]** :
  Supervised patch detection

### CellGraphAssay(5), PNAAssay(5) and CellGraph methods

- [`RenameCells(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md)
  [`show(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md)
  [`subset(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md)
  [`merge(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md) :
  CellGraphAssay Methods

- [`RenameCells(`*`<CellGraphAssay5>`*`)`](CellGraphAssay5-methods.md)
  [`show(`*`<CellGraphAssay5>`*`)`](CellGraphAssay5-methods.md)
  [`subset(`*`<CellGraphAssay5>`*`)`](CellGraphAssay5-methods.md)
  [`merge(`*`<CellGraphAssay5>`*`)`](CellGraphAssay5-methods.md)
  [`JoinLayers(`*`<CellGraphAssay5>`*`)`](CellGraphAssay5-methods.md) :
  CellGraphAssay5 Methods

- [`RenameCells(`*`<PNAAssay>`*`)`](PNAAssay-methods.md)
  [`show(`*`<PNAAssay>`*`)`](PNAAssay-methods.md)
  [`subset(`*`<PNAAssay>`*`)`](PNAAssay-methods.md)
  [`merge(`*`<PNAAssay>`*`)`](PNAAssay-methods.md) : PNAAssay Methods

- [`RenameCells(`*`<PNAAssay5>`*`)`](PNAAssay5-methods.md)
  [`show(`*`<PNAAssay5>`*`)`](PNAAssay5-methods.md)
  [`subset(`*`<PNAAssay5>`*`)`](PNAAssay5-methods.md)
  [`merge(`*`<PNAAssay5>`*`)`](PNAAssay5-methods.md)
  [`JoinLayers(`*`<PNAAssay5>`*`)`](PNAAssay5-methods.md) : PNAAssay5
  Methods

- [`RenameCells(`*`<MPXAssay>`*`)`](MPXAssay-methods.md)
  [`show(`*`<MPXAssay>`*`)`](MPXAssay-methods.md)
  [`subset(`*`<MPXAssay>`*`)`](MPXAssay-methods.md)
  [`merge(`*`<MPXAssay>`*`)`](MPXAssay-methods.md) : MPXAssay Methods

- [`CreateCellGraphAssay()`](CreateCellGraphAssay.md) : Create a
  CellGraphAssay object

- [`CreateCellGraphAssay5()`](CreateCellGraphAssay5.md) : Create a
  CellGraphAssay5 object

- [`CreatePNAAssay()`](CreatePNAAssay.md) : Create a PNAAssay object

- [`CreatePNAAssay5()`](CreatePNAAssay5.md) : Create a PNAAssay5 object

- [`CellGraphs()`](CellGraphs.md)
  [`` `CellGraphs<-`() ``](CellGraphs.md) : CellGraphs

- [`PolarizationScores()`](PolarizationScores.md)
  [`` `PolarizationScores<-`() ``](PolarizationScores.md) :
  PolarizationScores

- [`ColocalizationScores()`](ColocalizationScores.md)
  [`` `ColocalizationScores<-`() ``](ColocalizationScores.md) :
  ColocalizationScores

- [`ProximityScores()`](ProximityScores.md)
  [`` `ProximityScores<-`() ``](ProximityScores.md) : ProximityScores

- [`Edgelists()`](Edgelists.md) : Load edgelists

- [`Normalize()`](Normalize.md) : Normalize MPX or PNA data

- [`FSMap()`](FSMap.md) [`` `FSMap<-`() ``](FSMap.md) : FS map

- [`RestorePaths()`](RestorePaths.md) : Restore PXL file paths

- [`CellGraphData()`](CellGraphData.md)
  [`` `CellGraphData<-`() ``](CellGraphData.md) : Get and set CellGraph
  object data

- [`as.CellGraphAssay()`](as.CellGraphAssay.md) :

  Convert objects to a `CellGraphAssay`

- [`as.CellGraphAssay5()`](as.CellGraphAssay5.md) :

  Convert objects to a `CellGraphAssay5`

- [`as.PNAAssay()`](as.PNAAssay.md) :

  Convert objects to a `PNAAssay`

- [`as.PNAAssay5()`](as.PNAAssay5.md) :

  Convert objects to a `PNAAssay5`

- [`CellGraph-class`](CellGraph-class.md)
  [`CellGraph`](CellGraph-class.md) : The CellGraph class

- [`CreateCellGraphObject()`](CreateCellGraphObject.md) : Create a
  CellGraph object

### Graph functions

- [`edgelist_to_simple_Anode_graph()`](graph-conversion.md)
  **\[deprecated\]** : A-node projection

- [`edgelist_to_simple_bipart_graph()`](edgelist_to_simple_bipart_graph.md)
  : Create a simple bipartite graph from an edgelist

- [`node_markers_counts()`](node_markers_counts.md) **\[deprecated\]** :
  Calculate antibody counts per A-node

- [`color_by_marker()`](color_by_marker.md) **\[deprecated\]** :

  Add node colors to a `CellGraph`

- [`compute_transition_probabilities()`](compute_transition_probabilities.md)
  : Compute transition probabilities

- [`cos_distance_weights()`](edge-weights-pmds.md)
  [`prob_distance_weights()`](edge-weights-pmds.md) **\[experimental\]**
  : Calculate edge weights for pMDS

- [`local_G()`](local_G.md) : Calculate Local G

- [`expand_adjacency_matrix()`](expand_adjacency_matrix.md) : Expand an
  adjacency matrix to include higher-order neighborhoods

### Export

- [`WriteMPX_pxl_file()`](WriteMPX_pxl_file.md) **\[experimental\]** :
  Export Seurat object MPX data to a PXL file
- [`export_plot()`](export_plot.md) : Save a ggplot2 plot to one or more
  file formats with options for size and directory creation.

### Classes

- [`CellGraphAssay-class`](CellGraphAssay-class.md)
  [`CellGraphAssay`](CellGraphAssay-class.md) : The CellGraphAssay class
- [`CellGraphAssay5-class`](CellGraphAssay5-class.md)
  [`CellGraphAssay5`](CellGraphAssay5-class.md) : The CellGraphAssay5
  class
- [`MPXAssay-class`](MPXAssay-class.md) [`MPXAssay`](MPXAssay-class.md)
  : MPXAssay class
- [`PNAAssay-class`](PNAAssay-class.md) [`PNAAssay`](PNAAssay-class.md)
  : The PNAAssay class
- [`PNAAssay5-class`](PNAAssay5-class.md)
  [`PNAAssay5`](PNAAssay5-class.md) : The PNAAssay5 class

### Misc

- [`pixelatorR_options`](pixelatorR_options.md) : Global Options in
  pixelatorR
- [`pixelatorR-package`](pixelatorR-package.md) : pixelatorR: Data
  Structures, Data Processing Tools and Visualization Tools for MPX
  Single Cell Data
- [`minimal_mpx_pxl_file()`](mpx_dataset.md) : Five Cells MPX Test data
- [`minimal_pna_pxl_file()`](pna_dataset.md) : Five Cells PNA Test data
- [`show(`*`<CellGraph>`*`)`](CellGraph-methods.md)
  [`subset(`*`<CellGraph>`*`)`](CellGraph-methods.md) : CellGraph
  Methods
- [`RenameCells(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md)
  [`show(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md)
  [`subset(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md)
  [`merge(`*`<CellGraphAssay>`*`)`](CellGraphAssay-methods.md) :
  CellGraphAssay Methods
- [`inspect_pxl_file()`](inspect_pxl_file.md) : Inspect a PXL file with
  MPX data
- [`pixelatorR_style()`](pixelatorR_style.md) : The pixelatorR style
- [`print(`*`<pixelator_metadata>`*`)`](print.pixelator_metadata.md) :
  Print method for pixelator_metadata
- [`assert_single_value()`](type_check_helpers.md)
  [`assert_vector()`](type_check_helpers.md)
  [`assert_class()`](type_check_helpers.md)
  [`assert_mpx_assay()`](type_check_helpers.md)
  [`assert_pna_assay()`](type_check_helpers.md)
  [`assert_pixel_assay()`](type_check_helpers.md)
  [`assert_within_limits()`](type_check_helpers.md)
  [`assert_function()`](type_check_helpers.md)
  [`assert_file_exists()`](type_check_helpers.md)
  [`assert_file_ext()`](type_check_helpers.md)
  [`assert_x_in_y()`](type_check_helpers.md)
  [`assert_single_values_are_different()`](type_check_helpers.md)
  [`assert_singles_match()`](type_check_helpers.md)
  [`assert_length()`](type_check_helpers.md)
  [`assert_max_length()`](type_check_helpers.md)
  [`assert_vectors_x_y_length_equal()`](type_check_helpers.md)
  [`assert_unique()`](type_check_helpers.md)
  [`assert_col_class()`](type_check_helpers.md)
  [`assert_col_in_data()`](type_check_helpers.md)
  [`assert_non_empty_object()`](type_check_helpers.md)
  [`assert_is_one_of()`](type_check_helpers.md)
  [`assert_different()`](type_check_helpers.md)
  [`assert_vectors_match()`](type_check_helpers.md)
  [`assert_valid_color()`](type_check_helpers.md) : Type check helpers
- [`isotype_pls()`](isotype_pls.md) : Partial Least Squares Regression
  for background factor regression
- [`get_duckdb_config()`](get_duckdb_config.md) : Fetches the config
  values to be used with DuckDB
- [`CalculateDispersion()`](CalculateDispersion.md) : Calculate
  dispersion of components
