# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2024-??-??

### Updates 

- The R arrow version is no longer pinned to v14. This allows the package to be installed with the latest version of arrow.
- Updated `LoadCellGraphs` methods to be compatible with R arrow v17
- `RunDPA` and `RunDCA` now handles multiple `targets` for differential tests. Previously, only 1 `target` could be compared against `reference`. Now, if multiple `targets` are provided, the function will perform multiple differential tests, one for each `target` against `reference`. This is typically useful when comparing multiple conditions against a single control group.
- `ColocalizationHeatmap` has been made more flexible, such that any column names in the input data can be used as long as the data has a data format suitable for a heat map. 

### Fixes

- In `ColocalizationScoresToAssay`: Changed marker pair separator from "-" to "/", to avoid string operation issues due to "-" occurring in marker names.

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
