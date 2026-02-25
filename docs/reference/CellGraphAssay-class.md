# The CellGraphAssay class

The CellGraphAssay object is an extended
[`Assay`](https://satijalab.github.io/seurat-object/reference/Assay-class.html)
for the storage and analysis of MPX single-cell data.

## Slots

- `cellgraphs`:

  A named list of [`CellGraph`](CellGraph-class.md) objects

- `polarization`:

  A `tbl_df` with polarization scores

- `colocalization`:

  A `tbl_df` with colocalization scores

- `fs_map`:

  A `tbl_df` with information on source pxl file paths, sample IDs, and
  component IDs
