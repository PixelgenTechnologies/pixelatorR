# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added
- Supervised patch detection implemented in the `patch_detection` function. Added `identify_markers_for_patch_analysis` to identify markers for patch analysis.
- `render_rotating_layout` function to create videos of rotating cells from a tibble containing layout coordinates. The function supports multiple video formats but uses GIF as default.
- `subset` method for `CellGraph` class
- Option to add marker count proportions to the proximity score table in `ProximityScores`
- `pack_2bits` and `unpack_2bits` to pack and unpack DNA sequences into 64-bit integers using 2 bits per base.
- Experimental `PredictDoublets` function for detecting doublets in a Seurat object or count matrix.
- `SimulateDoublets` to simulate doublets.
- `FindAnnoyNeighbors` Computes nearest neighbors using the Annoy algorithm.
- `DensityScatterPlot` now has an argument `equal_axes` to control whether the x and y axes should have a common range. 
- `return_id` argument to `SimulateDoublets` to output the IDs of cells used to simulate each doublet.
- `PredictDoublets` can now be run iteratively using the `iter` argument to increase robustness.

### Updated
- `PredictDoublets` now have `simulation_rate = 1` and  `p_threshold = 0.05` as defaults.

### Fixes
- Fixed bug in `ColocalizationHeatmap` where `marker1_col` and `marker2_col` only worked for "marker_1" and "marker_2".
- Fixed bug in `DensityScatterPlot` where the % cells label would be calculated across all facets instead of per each facet.

## [0.13.0] - 2025-04-23

### Added
- `PixelDB` R6 class to access data from a PXL file (`<pxl_file>`) containing a duckdb database.
  - `PixelDB$new` create a new `PixelDB` object from a PXL file containing PNA data.
  - `PixelDB$info` get information about the tables stored in the PXL file.
  - `PixelDB$query` send an SQL query to the database.
  - `PixelDB$check_connection` check if the connection to the PXL file is still valid.
  - `PixelDB$reconnect` reconnect to the database if the connection is closed.
  - `PixelDB$names` get the names of the tables stored in the database.
  - `PixelDB$fetch_table` fetch an entire table as a `data.frame`.
  - `PixelDB$fetch_table_subset` fetch a subset of a table as a `data.frame`
  - `PixelDB$counts` fetch the antibody count matrix.
  - `PixelDB$proximity` fetch the proximity scores table.
  - `PixelDB$cell_meta` fetch the component/cell meta data.
  - `PixelDB$protein_meta` fetch the protein meta data.
  - `PixelDB$run_meta` fetch the Pixelator run meta data.
  - `PixelDB$components_edgelist` fetch the edgelist(s) for selected components/cells.
  - `PixelDB$components_layout` fetch the layout(s) for selected components/cells.
  - `PixelDB$components_marker_counts` fetch the node counts for selected component/cell graphs.
  - `PixelDB$export_parquet` export a table in the database to a parquet file.
  - `PixelDB$close` close the connection.
- `ReadPNA_counts` function to load the count matrix from a PXL file containing PNA data.
- `ReadPNA_proximity` function to load the proximity scores table from a PXL file containing PNA data. Also supports lazy loading.
- `ReadPNA_edgelist` function to load the edgelist from a PXL file containing PNA data. Also supports lazy loading.
- `ReadPNA_layouts` function to load component layouts from a PXL file containing PNA data with pre-computed layouts.
- `ReadPNA_Seurat` function to construct a `Seurat` object from a PXL file containing PNA data.
- `ReadPNA_metadata` function to load sample meta data from a PXL file containing PNA data.
- `PNAAssay` class to store PNA data in a `Seurat` object (v3).
- `CreatePNAAssay` to create a `PNAAssay` object.
- `PNAAssay5` class to store PNA data in a `Seurat` object (v5).
- `CreatePNAAssay5` to create a `PNAAssay5` object.
- `Edgelists` methods for `PNAAssay`, `PNAAssay5` and `Seurat` to load edgelists. Supports lazy loading for manipulation with `dbplyr`.
- `ProximityScores` methods for `PNAAssay`, `PNAAssay5` and `Seurat` to fetch proximity scores. Supports lazy loading for manipulation with `dbplyr`.
- `ProximityScores<-` methods for `PNAAssay`, `PNAAssay5` and `Seurat` to set proximity scores.
- `ProximityScoresToAssay` methods for `data.frame`, `tbl_lazy`, `PNAAssay`, `PNAAssay5` and `Seurat` to convert the long formatted proximity score table into a wide format.
- `LoadCellGraphs` methods for for `PNAAssay` and `PNAAssay5`.
- `ComputeLayout` methods for for `PNAAssay` and `PNAAssay5`.
- `RemoveCellGraphs` methods for for `PNAAssay` and `PNAAssay5`.
- `CellGraphs` methods for for `PNAAssay` and `PNAAssay5`.
- `RestorePaths` methods for for `PNAAssay` and `PNAAssay5`.
- `FSMap`/`FSMap<-` methods for `PNAAssay` and `PNAAssay5`.
- `show` method for `PNAAssay` and `PNAAssay5`.
- `subset` method for `PNAAssay` and `PNAAssay5`.
- `merge` method for `PNAAssay` and `PNAAssay5`.
- `RenameCells` methods for `PNAAssay` and `PNAAssay5`.
- `JoinLayers` method for `PNAAssay5`.
- `as.PNAAssay` method to convert an `Assay` object to a `PNAAssay` object.
- `as.PNAAssay5` method to convert an `Assay5` object to a `PNAAssay5` object.
- `DifferentialProximityAnalysis` function to perform differential testing on PNA proximity scores. The function has a similar API as `RunDAA`, `RunDPA` and `RunDCA` but uses a much faster implementation of the Wilcoxon rank sum test (Mann-Whitney U test) with the `data.table` R package.
- A minimal PXL file woth PNA data.
- `minimal_mpx_pxl_file` function to get the path to the minimal MPX PXL file.
- `minimal_mpx_pna_file` function to get the path to the minimal PNA PXL file.
- Utility function for asserting valid colors; `assert_valid_color`.
- Utility function for asserting valid `PNAAssay`/`PNAAssay5`; `assert_pna_assay`.
- Utility function for asserting valid `PNAAssay`/`PNAAssay5`/`CellGraphAssay`/`CellGraphAssay5`; `assert_pna_assay`.

### Updates
- `MoleculeRankPlot` now supports Seurat objects with PNA data with `n_umi` representing the total number of detected antibodies.
- `TauPlot` now supports Seurat objects with PNA data using `n_umi` on the y-axis.
- `Plot2DGraph` now supports Seurat objects with PNA data. 
- `Plot2DGraphM` now supports Seurat objects with PNA data. 
- `Plot3DGraph` now supports Seurat objects with PNA data. 
- `DensityScatterPlot` can now draw `rectangle` or `quadrant` gates by selecting the appropriate `gate_type` argument. Additionally, gate annotation aesthetics can now be customized using `annotation_params`.

### Changes
- `ComputeLayout` now only supports the "pmds" and "wpmds" graph drawing methods. The "kk", "fr" and "drl" methods have been removed but can be run if needed using the `custom_layout_function` parameter. The default layout method is now "wpmds" with `dim = 3`.
- `NormalizeMPX` is superseded by `Normalize`. The `NormalizeMPX` function will be removed in a future release.

### Removed
- `LoadCellGraphs.data.frame` method

### Fixes
- Fixed bug in `DensityScatterPlot` where the `gate_type` default would lead to an error.
- Fixed bug in `DensityScatterPlot` where the x- and y-axis titles were hardcoded as "Marker1" and "Marker2"
- Fixed bug in `subset.MPXAssay` and `subset.PNAAssay` where the `fs_map` table was not filtered correctly when all components from a sample are removed.

## [0.12.1] - 2025-01-21

### Fixes
- `MoleculeRankPlot` now also accepts columns "edges" and "molecules" as numeric class, in addition to integer class.
- `ColocalizationHeatmap` now has the same order of x and y axis when `symmetrise = TRUE` and `type = "dots"`.
- Updated ".lintr" rules to handle return rule added in lintr 3.2.0.
- Updated deprecated v3 of GitHub Action `upload-artifact` to v4.
- Swapped `zip::unzip` with `utils::unzip` in `ReadMPX_counts`

## [0.12.0] - 2025-01-16

### Added

- `RunDAA` : Differential abundance analysis function with a similar interface to `RunDPA` and `RunDCA`. `RunDAA` uses the `FindMarkers` function from Seurat to perform differential abundance analysis, but enables splitting of tests into multiple groups. By default, it reports the difference in means instead of `avg_log2FC`.

### Updates

- Updated type assertions and improved error messaging (inspired by the tidyverse style guide).
- The `ComputeLayout.Seurat` method now supports parallelized computation of layouts.
- Added option to fetch marker counts in `PolarizationScores` and `ColocalizationScores` methods. This is for example useful when filtering spatial metrics tables for markers with low counts.
- Silenced warnings in `RunDPA`/`RunDCA` when running tests in parallel to avoid halting the R session.
- Improved clean up of temporary files created by `ReadMPX_counts` and `ReadMPX_item`. 
- `RunDPA` and `RunDCA` now accepts any numeric vector from the spatial metric table as input for differential testing. The metric is specified by `polarity_metric` (`RunDPA`) or `coloc_metric` (`RunDCA`).
- Updated `subset` and `merge` methods for `MPXAssay` to have less stringent validation of the spatial metric tables (polarity and colocalization scores).  
- Updated `ReadMPX_Seurat` to have less stringent validation of the spatial metric tables (polarity and colocalization scores).
- `ColocalizationHeatmap` now allows legend titles and the legend range to be manually set

### Fixes

- pixelatorR read functions now uses `utils::unzip` instead of `zip::unzip` to support PXL files larger than 2GB
- `LoadCellGraphs` now throws an error if duplicated cell ids (`cells`) are provided
- PXL files missing spatial scores can now be loaded with `ReadMPX_Seurat` without throwing an error. This is useful when the pixelator pipeline was run without computing spatial scores.

## [0.11.0] - 2024-09-18

### Updates 

- The R arrow version is no longer pinned to v14. This allows the package to be installed with the latest version of arrow.
- Updated `LoadCellGraphs` methods to be compatible with R arrow v17
- `RunDPA` and `RunDCA` now handles multiple `targets` for differential tests. Previously, only 1 `target` could be compared against `reference`. Now, if multiple `targets` are provided, the function will perform multiple differential tests, one for each `target` against `reference`. This is typically useful when comparing multiple conditions against a single control group.
- `ColocalizationHeatmap` has been made more flexible, such that any column names in the input data can be used as long as the data has a data format suitable for a heat map. 

### Fixes

- In `ColocalizationScoresToAssay`: Changed marker pair separator from "-" to "/", to avoid string operation issues due to "-" occurring in marker names.
- In `DensityScatterPlot`: Fixed a bug where the reported percentage of cells in a gate was incorrect when using multiple gates. 
- `NormalizeMPX` now also handles Assay and Assay5 objects.

### Added

- `RestorePaths` : updates the PXL file paths in a `CellGraphAssay`, a `CellGraphAssay5` or a `Seurat` object created with pixelatorR. This function is useful when PXL files have been moved to a different location or when sharing Seurat objects with other users which would cause `LoadCellGraphs` to fail.
- `ReadMPX_metadata` : loads metadata from a PXL file. A `print` method for the output returned by `ReadMPX_metadata` is also included to provide a summary of the metadata.

### [0.10.2] - 2024-07-24

### Fixes

- fixed bug in `TauPlot` where `group_by` was not properly evaluated

### [0.10.1] - 2024-07-23

### Fixes

- updated `TauPlot` to handle new metric names introduced in pixelator v0.18. For data produced with pixelator v0.18, `TauPlot` now uses `mean_molecules_per_a_pixel` instead of `umi_per_upia`.
- Updated `DensityScatterPlot` to not use deprecated `dplyr` functionality. 
- Changed `ComputeLayout` to add a suffix ('_3d') to the layout name when `dim = 3`. This makes the naming of layouts consistent with the naming used in pixelator (Python).

## [0.10.0] - 2024-07-10

### Added

- `layout_with_weighted_pmds` : Compute a graph layout using weighted PMDS.
- option to use `layout_with_weighted_pmds` in `ComputeLayout` by setting `layout_method = "wpmds"`

### Fix

- Updated documentation for `local_G`. The equations for the Z-scores have been corrected according to the original publication by Ord and Getis.

### Added 

- `AbundanceColocalizationPlot` : Plot a scatter plot of the abundance of two sets of markers, colored by the 
colocalization score the marker pairs have in each component. 

## [0.9.0] - 2024-07-03

### Added 

- `NormalizeMPX` : Normalize MPX data using either dsb or CLR transformations.

### Updates 

- `PseudoDensityPlot` has been renamed to `DensityScatterPlot` to better reflect the function's purpose and 
documentation has been improved

## [0.8.1] - 2024-07-01

### Fixes

- Fixed bug where `load_polarity_scores` controls whether colocalization scores are loaded instead of `load_colocalization_scores` in `ReadMPX_Seurat`.

## [0.8.0] - 2024-06-28

### Added

- `pseudocolor_plot` plotting function to create scatterplots colored by density.

## [0.7.2] - 2024-06-27

### Fixes

- Updated `ReadMPX_layouts` to handle new layout parquet file format (including node names). This update is necessary to ensure that the order of pre-computed layout node coordinates match the order of nodes in loaded graphs.
- Updated `LoadCellGraphs` to add pre-computed layouts to `CellGraph` objects with correct node order.
- Updated `WriteMPX_pxl_file` to include node names in the exported layout parquet files. 

## [0.7.1] - 2024-06-24

### Added

- `FSMap` Seurat method

## [0.7.0] - 2024-06-13

### Added

- `MoleculeRankPlot` : replaces `EdgeRankPlot`. `EdgeRankPlot` is now deprecated and will be removed in a future release.

## [0.6.1] - 2024-06-11

### Fixes

- Fixed bug in `local_G` when `use_weights=FALSE` and `k>1`

## [0.6.0] - 2024-06-10

### Added

- `local_G` : compute Z-scores for marker counts in a component graph comparing local vs global abundance
- `JoinLayers.CellGraphAssay5` : `JoinLayers` method for `CellGraphAssay5` class objects to support joining layers in Seurat V5 type Assay objects supported by pixelatorR

## [0.5.2] - 2024-05-20

### Fixes

- Fixed bug in `.validate_polarization` that would throw an error if columns were in unexpected order in `polarization.parquet`. 

## [0.5.1] - 2024-05-20

### Updates

- Loading A-node projected graphs from PXL files with upia/upib stored as dictionaries is now supported when using `LoadCellGraphs`

## [0.5.0] - 2024-05-16

### Added

- `inspect_pxl_file` : utility function to inspect the contents of a PXL file
- `ReadMPX_layouts` : read function to import pre-computed layouts

### Updates

- Updated show method for `CellGraph` class objects

### Fixes

- Updated experimental function `WriteMPX_pxl_file` for compatibility with the anndata Python library. Note that since the `WriteMPX_pxl_file` function is experimental, the function may change at any point in the future.

## [0.4.2] - 2024-05-06

### Fixes

- Fixed bug in `PolarizationSsoresToAssay` and `PolarizationSsoresToAssay` which failed when column IDs included dashes

## [0.4.1] - 2024-04-29

### Updates

- Seurat methods `PolarizationSsoresToAssay` and `PolarizationSsoresToAssay` now adds a key to the returned assay object

## [0.4.0] - 2024-04-24

### Updated `pixelatorR` classes and their methods for improved I/O

- `pixelatorR` now provides two classes to store MPX data: `CellGraphAssay5` and `CellGraphAssay5`
  - `CellGraphAssay5` inherits the `Assay5` class introduced in Seurat v5 and will be used when `options(Seurat.object.assay.version = "v5")`
  - `CellGraphAssay` inherits the `Assay` class from Seurat v3 and will be used when `options(Seurat.object.assay.version = "v3")`
- All `CellGraphAssay` methods now handle `CellGraphAssay5` class objects
  
### hdf5

- Switched from R package `rhdf5` (Bioconductor) to `hdf5r` (CRAN) to handle reading and writing of HDF5 files

### Writing PXL files

- Introduced experimental function `WriteMPX_pxl_file` to write MPX data from a Seurat object to a PXL file

### Fixes

- Made `EdgeRankPlot` compatible with incoming changes in pixelator Python >0.16.2 where `edges` is renamed to `molecules` in `CellGraphAssay` objects

## [0.3.0] 2024-03-28

### Added utility functions to clean up edgelist directories

- `edgelist_directories_clean` runs a cleanup of the edgelist directory set by the "pixelatorR.arrow_outdir" global option
- `edgelist_directories_du` returns the total disk usage of the edgelist directory set by the "pixelatorR.arrow_outdir" global option
- Added global option "pixelatorR.arrowdir_maxsize". When running `ReadMPX_Seurat`, `merge.CellGraphAssay`, `subset.CellGraphAssay` or `RenameCells.CellGraphAssay`, a clean up will be triggered if the total disk usage of the edgelist directory exceeds this value. The default value is 5 GB.

## [0.2.0] - 2024-03-18

### Added new test data set - five_cells.pxl

Once `pixelatorR` is installed, the test data can be read from the following path:

````
pxl_file <- system.file("extdata/five_cells",
  "five_cells.pxl",
  package = "pixelatorR"
)
````

- The data set contains 5 cells from three different MPX data sets:
   - [CD3 capped cell: RCVCMP0000217](https://software.pixelgen.com/datasets/) from the CD3 capping data set (stimulated)
   - [B cell: RCVCMP0000118, CD4 T cell: RCVCMP0000487, CD8 T cell: RCVCMP0000655](https://software.pixelgen.com/datasets/) from a PBMC data set
   - [Uropod cell: RCVCMP0000263](https://software.pixelgen.com/datasets/) from the RANTES stimulated data set

- All cells now have fully connected components
- The total file size is 2.1 MB, where the following fields have been excluded from the edgelist.parquet file: "umi", "sequence", "umi_unique_count", "upi_unique_count"

## [0.1.3] - 2024-03-18

- The latest version of `arrow` on CRAN (v15.0.1) contains several bugs which are picked up by the R CMD checks. Fixed the version of `arrow` to v14.0.2.1 until a stable release is available. 

## [0.1.2] - 2024-03-14

* Fixed bug in `LoadCellGraphs` when working on Seurat objects with multiple data sets (merged). Previous version returned components in an invalid order.

## [0.1.1] - 2024-02-29

* Fixed bug in `ComputeLayout` when normalizing or projecting layout coordinates
* Added layout utility functions:
  * `center_layout_coordinates` - Center layout coordinates at the origin
  * `normalize_layout_coordinates` - Normalize layout coordinates such that they are centered at origin and the have a median length (euclidean norm) of 1
  * `project_layout_coordinates_on_unit_sphere` - Project layout 3D coordinates onto the unit sphere

## [0.1.0] - 2024-02-21

* First public release of pixelatorR.
