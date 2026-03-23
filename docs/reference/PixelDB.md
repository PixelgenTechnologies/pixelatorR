# PXL database class

PXL database class

PXL database class

## Details

This class provides an interface for working with a PXL file. A PXL file
is a duckdb database that contains various tables created by the
Pixelator data processing pipeline.

The class provides methods to query the database and extract information
from it. The class is intended to be used internally by the package but
is exposed to the user for advanced use cases.

Example usage:

    # Create a PixelDB object with a connection to a PXL file
    pxl_db <- PixelDB$new("path/to/pxl/file")
    pxl_db$info()

    # Fetch count matrix
    counts <- pxl_db$counts()

    # Fetch proximity scores
    proximity <- pxl_db$proximity()

## Methods

### Public methods

- [`PixelDB$new()`](#method-PixelDB-new)

- [`PixelDB$info()`](#method-PixelDB-info)

- [`PixelDB$query()`](#method-PixelDB-query)

- [`PixelDB$reconnect()`](#method-PixelDB-reconnect)

- [`PixelDB$check_connection()`](#method-PixelDB-check_connection)

- [`PixelDB$names()`](#method-PixelDB-names)

- [`PixelDB$fetch_table()`](#method-PixelDB-fetch_table)

- [`PixelDB$fetch_table_subset()`](#method-PixelDB-fetch_table_subset)

- [`PixelDB$counts()`](#method-PixelDB-counts)

- [`PixelDB$proximity()`](#method-PixelDB-proximity)

- [`PixelDB$cell_meta()`](#method-PixelDB-cell_meta)

- [`PixelDB$protein_meta()`](#method-PixelDB-protein_meta)

- [`PixelDB$run_meta()`](#method-PixelDB-run_meta)

- [`PixelDB$components_edgelist()`](#method-PixelDB-components_edgelist)

- [`PixelDB$components_layout()`](#method-PixelDB-components_layout)

- [`PixelDB$components_marker_counts()`](#method-PixelDB-components_marker_counts)

- [`PixelDB$export_parquet()`](#method-PixelDB-export_parquet)

- [`PixelDB$close()`](#method-PixelDB-close)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Set up a connection to a PXL file

#### Usage

    PixelDB$new(file)

#### Arguments

- `file`:

  A path to a PXL file

#### Returns

An `R6` object representing the PXL database

#### Examples

    library(dplyr)

    pxl_file <- minimal_pna_pxl_file()
    db <- PixelDB$new(pxl_file)

------------------------------------------------------------------------

### Method `info()`

Show information about tables in the PXL file

#### Usage

    PixelDB$info()

#### Returns

A `tbl_df` with information about tables in the PXL file

#### Examples

    db <- PixelDB$new(pxl_file)
    db$info()

------------------------------------------------------------------------

### Method `query()`

Send a database query

#### Usage

    PixelDB$query(sql)

#### Arguments

- `sql`:

  An SQL query

#### Returns

A `data.frame` with the results of the query

#### Examples

    # Select the proximity score table
    db$query("SELECT * FROM proximity") %>% head()

------------------------------------------------------------------------

### Method `reconnect()`

Method to reconnect to the PXL file. This is useful if the connection is
lost, for instance if the R session restarted.

If the method fails, an error is thrown.

#### Usage

    PixelDB$reconnect()

#### Returns

Nothing

#### Examples

    db$close()
    db$reconnect()

------------------------------------------------------------------------

### Method `check_connection()`

Check the connection to the PXL database

If the connection is invalid, the `$reconnect()` method is called.

#### Usage

    PixelDB$check_connection()

#### Returns

Nothing

#### Examples

    # If the connection is closed, the method will attempt to reconnect
    db$check_connection()

------------------------------------------------------------------------

### Method [`names()`](https://rdrr.io/r/base/names.html)

Get the names of the tables in the database

#### Usage

    PixelDB$names()

#### Returns

A character vector with table names that can be used in `$fetch_table()`

#### Examples

    # Get the table names
    db$names()

------------------------------------------------------------------------

### Method `fetch_table()`

Fetches an entire table from the database

This general method can be used to fetch any table from the database.
Large tables such as the edgelist are not recommended to be fetched in
this way as you will get the entire table in memory.

#### Usage

    PixelDB$fetch_table(name)

#### Arguments

- `name`:

  The name of the table

#### Returns

A `data.frame` with the table contents

#### Examples

    # Fetch any table from the database
    db$fetch_table("proximity") %>% head()

------------------------------------------------------------------------

### Method `fetch_table_subset()`

Fetches a subset of columns from a table in the database

This method can be used to fetch a subset of columns from a table in the
database.

#### Usage

    PixelDB$fetch_table_subset(name, columns_filter)

#### Arguments

- `name`:

  The name of the table

- `columns_filter`:

  A named list

#### Returns

A `data.frame` with the table contents

#### Examples

    # Fetch any table from the database and filter on the fly
    db$fetch_table_subset(
      "proximity",
      columns_filter = list("component" = c("0a45497c6bfbfb22", "2708240b908e2eba"))
    ) %>% head()

------------------------------------------------------------------------

### Method `counts()`

Fetches the \_\_adata\_\_X count table from the database and converts it
to a sparse `dgCMatrix`

#### Usage

    PixelDB$counts()

#### Returns

A sparse matrix of class `dgCMatrix` with antibody counts

#### Examples

    # Fetch the antibody counts
    X <- db$counts()
    X[1:4, 1:4]

------------------------------------------------------------------------

### Method `proximity()`

Fetches the proximity scores from the database

#### Usage

    PixelDB$proximity(calc_log2_ratio = TRUE, lazy = FALSE)

#### Arguments

- `calc_log2_ratio`:

  A logical specifying whether to calculate and add a log2ratio column
  to the output table. Default is `TRUE`

- `lazy`:

  A logical specifying whether to load the data lazily. If `TRUE`, a
  `tbl_lazy` object is returned.

#### Returns

A `tbl_df` with the proximity scores

#### Examples

    # Fetch the proximity scores
    prox <- db$proximity()
    prox %>% head()

------------------------------------------------------------------------

### Method `cell_meta()`

Fetches the \_\_adata\_\_obs meta data

#### Usage

    PixelDB$cell_meta()

#### Returns

A `data.frame` with the cell meta data

#### Examples

    # Fetch the cell meta data
    db$cell_meta() %>% head()

------------------------------------------------------------------------

### Method `protein_meta()`

Fetches the \_\_adata\_\_var meta data

#### Usage

    PixelDB$protein_meta()

#### Returns

A `data.frame` with the protein meta data

#### Examples

    # Fetch the protein meta data
    db$protein_meta() %>% head()

------------------------------------------------------------------------

### Method `run_meta()`

Fetch the run meta data

#### Usage

    PixelDB$run_meta()

#### Returns

A `tbl_df` with the run meta data

#### Examples

    # Fetch the run meta data
    db$run_meta()

------------------------------------------------------------------------

### Method `components_edgelist()`

Fetches a component edgelist or the entire edgelist from the database

The UMIs are encoded as int64 but since R doesn't support int64, the
UMIs can be converted to character vectors be specifying the
`umi_data_type`.

#### Usage

    PixelDB$components_edgelist(
      components,
      umi_data_type = c("int64", "string", "suffixed_string"),
      lazy = FALSE,
      include_all_columns = FALSE
    )

#### Arguments

- `components`:

  A character vector with component names or NULL to get all components

- `umi_data_type`:

  One of "int64", "string" or "suffixed_string". Default is "int64".

  - "int64": The UMIs are encoded as int64

  - "string": The UMIs are encoded as character

  - "suffixed_string": The UMIs are encoded as character with a suffix
    '-umi1' or '-umi2' added

- `lazy`:

  A logical specifying whether to load the data lazily. If `TRUE`, a
  `tbl_lazy` object is returned.

- `include_all_columns`:

  Logical specifying whether to include all columns in the output.

#### Returns

A `data.frame` with the component edgelist:

- umi1: A unique ID of the first RCA product

- umi2: A unique ID of the second RCA product

- marker_1: The first protein

- marker_2: The second protein

- component: The component name

- read_count: The number of reads supporting the edge

- uei_count: The number of unique event identifiers (UEIs) supporting
  the edge

#### Examples

    # Fetch edgelists
    db$components_edgelist("0a45497c6bfbfb22") %>% head()

------------------------------------------------------------------------

### Method `components_layout()`

Fetches layout for selected components

This method fetches the x, y, z layout coordinates for selected PNA
`components`.

#### Usage

    PixelDB$components_layout(
      components,
      add_marker_counts = FALSE,
      verbose = TRUE
    )

#### Arguments

- `components`:

  A vector of component names or NULL for all components

- `add_marker_counts`:

  Add marker counts in wide format to the layout tables

- `verbose`:

  Print messages

#### Returns

A list with `tbl_df`'s with the layout coordinates and optionally marker
counts

#### Examples

    # Fetch layouts (NOTE: This will only work if the layouts exist in the database)
    \dontrun{
    db$components_layout("0a45497c6bfbfb22")[[1]] %>% head()
    }

------------------------------------------------------------------------

### Method `components_marker_counts()`

Fetches marker counts for components

This method fetches the node marker counts for selected `components`.
The node IDs are stored in the `name` column. If `components` is `NULL`,
all component marker counts are returned with an additional `components`
column in the resulting table.

#### Usage

    PixelDB$components_marker_counts(
      components,
      as_sparse = FALSE,
      verbose = FALSE
    )

#### Arguments

- `components`:

  A vector of component IDs or NULL for all components

- `as_sparse`:

  Return the marker counts as a sparse matrix

- `verbose`:

  Print messages

#### Returns

A list with `dgCMatrix` matrices or `tbl_df`'s with the marker counts

#### Examples

    # Fetch marker counts
    db$components_marker_counts("0a45497c6bfbfb22")[[1]][1:3, 1:4]

------------------------------------------------------------------------

### Method `export_parquet()`

Export a table to a parquet file

#### Usage

    PixelDB$export_parquet(
      parquet_file,
      table_name = c("proximity", "edgelist", "layouts"),
      compression = c("snappy", "zstd"),
      compression_level = 1L
    )

#### Arguments

- `parquet_file`:

  Path to the parquet file

- `table_name`:

  The name of the table to export

- `compression`:

  The compression algorithm to use. Options are 'snappy' and 'zstd'.

- `compression_level`:

  The compression level to use. Default is 1. Only used when compression
  is 'zstd'.

#### Returns

Nothing

#### Examples

    # Export a table to a parquet file
    tmp_parquet_file <- fs::file_temp(ext = "parquet")
    db$export_parquet(tmp_parquet_file, "proximity")
    fs::file_exists(tmp_parquet_file)

------------------------------------------------------------------------

### Method [`close()`](https://rdrr.io/r/base/connections.html)

Close connection

#### Usage

    PixelDB$close()

#### Returns

Nothing

#### Examples

    # Close the connection when finished
    db$close()

## Examples

``` r
## ------------------------------------------------
## Method `PixelDB$new`
## ------------------------------------------------

library(dplyr)

pxl_file <- minimal_pna_pxl_file()
db <- PixelDB$new(pxl_file)


## ------------------------------------------------
## Method `PixelDB$info`
## ------------------------------------------------

db <- PixelDB$new(pxl_file)
db$info()
#> # A tibble: 10 × 6
#>    database         schema name              column_names column_types temporary
#>    <chr>            <chr>  <chr>             <list>       <list>       <lgl>    
#>  1 minimal_PNA_PBMC main   __adata__X        <chr [159]>  <chr [159]>  FALSE    
#>  2 minimal_PNA_PBMC main   __adata__obs      <chr [24]>   <chr [24]>   FALSE    
#>  3 minimal_PNA_PBMC main   __adata__obsm_clr <chr [159]>  <chr [159]>  FALSE    
#>  4 minimal_PNA_PBMC main   __adata__obsm_lo… <chr [159]>  <chr [159]>  FALSE    
#>  5 minimal_PNA_PBMC main   __adata__uns      <chr [1]>    <chr [1]>    FALSE    
#>  6 minimal_PNA_PBMC main   __adata__var      <chr [6]>    <chr [6]>    FALSE    
#>  7 minimal_PNA_PBMC main   edgelist          <chr [7]>    <chr [7]>    FALSE    
#>  8 minimal_PNA_PBMC main   layouts           <chr [8]>    <chr [8]>    FALSE    
#>  9 minimal_PNA_PBMC main   metadata          <chr [1]>    <chr [1]>    FALSE    
#> 10 minimal_PNA_PBMC main   proximity         <chr [8]>    <chr [8]>    FALSE    


## ------------------------------------------------
## Method `PixelDB$query`
## ------------------------------------------------

# Select the proximity score table
db$query("SELECT * FROM proximity") %>% head()
#>   marker_1 marker_2 join_count join_count_expected_mean join_count_expected_sd
#> 1     CD56     CD56          0                     0.00              0.0000000
#> 2     CD56   mIgG2b          0                     0.03              0.1714466
#> 3     CD56     CD71          0                     0.01              0.1000000
#> 4     CD56      CD6          0                     2.07              1.5843617
#> 5     CD56 Siglec-9          0                     0.03              0.1714466
#> 6     CD56    CD79a          0                     0.00              0.0000000
#>   join_count_z join_count_p        component
#> 1      0.00000   0.50000000 c3c393e9a17c1981
#> 2     -0.03000   0.48803353 c3c393e9a17c1981
#> 3     -0.01000   0.49601064 c3c393e9a17c1981
#> 4     -1.30652   0.09568792 c3c393e9a17c1981
#> 5     -0.03000   0.48803353 c3c393e9a17c1981
#> 6      0.00000   0.50000000 c3c393e9a17c1981


## ------------------------------------------------
## Method `PixelDB$reconnect`
## ------------------------------------------------

db$close()
db$reconnect()
#> ✔ Connected


## ------------------------------------------------
## Method `PixelDB$check_connection`
## ------------------------------------------------

# If the connection is closed, the method will attempt to reconnect
db$check_connection()


## ------------------------------------------------
## Method `PixelDB$names`
## ------------------------------------------------

# Get the table names
db$names()
#>  [1] "__adata__X"          "__adata__obs"        "__adata__obsm_clr"  
#>  [4] "__adata__obsm_log1p" "__adata__uns"        "__adata__var"       
#>  [7] "edgelist"            "layouts"             "metadata"           
#> [10] "proximity"          


## ------------------------------------------------
## Method `PixelDB$fetch_table`
## ------------------------------------------------

# Fetch any table from the database
db$fetch_table("proximity") %>% head()
#>   marker_1 marker_2 join_count join_count_expected_mean join_count_expected_sd
#> 1     CD56     CD56          0                     0.00              0.0000000
#> 2     CD56   mIgG2b          0                     0.03              0.1714466
#> 3     CD56     CD71          0                     0.01              0.1000000
#> 4     CD56      CD6          0                     2.07              1.5843617
#> 5     CD56 Siglec-9          0                     0.03              0.1714466
#> 6     CD56    CD79a          0                     0.00              0.0000000
#>   join_count_z join_count_p        component
#> 1      0.00000   0.50000000 c3c393e9a17c1981
#> 2     -0.03000   0.48803353 c3c393e9a17c1981
#> 3     -0.01000   0.49601064 c3c393e9a17c1981
#> 4     -1.30652   0.09568792 c3c393e9a17c1981
#> 5     -0.03000   0.48803353 c3c393e9a17c1981
#> 6      0.00000   0.50000000 c3c393e9a17c1981


## ------------------------------------------------
## Method `PixelDB$fetch_table_subset`
## ------------------------------------------------

# Fetch any table from the database and filter on the fly
db$fetch_table_subset(
  "proximity",
  columns_filter = list("component" = c("0a45497c6bfbfb22", "2708240b908e2eba"))
) %>% head()
#>   marker_1 marker_2 join_count join_count_expected_mean join_count_expected_sd
#> 1    CD278    CD278          0                        0                      0
#> 2    CD278   HLA-DQ          0                        0                      0
#> 3    CD278    CD314          0                        0                      0
#> 4    CD278    VISTA          0                        0                      0
#> 5    CD278    CD326          0                        0                      0
#> 6    CD278     CD28          0                        0                      0
#>   join_count_z join_count_p        component
#> 1            0          0.5 0a45497c6bfbfb22
#> 2            0          0.5 0a45497c6bfbfb22
#> 3            0          0.5 0a45497c6bfbfb22
#> 4            0          0.5 0a45497c6bfbfb22
#> 5            0          0.5 0a45497c6bfbfb22
#> 6            0          0.5 0a45497c6bfbfb22


## ------------------------------------------------
## Method `PixelDB$counts`
## ------------------------------------------------

# Fetch the antibody counts
X <- db$counts()
X[1:4, 1:4]
#> 4 x 4 sparse Matrix of class "dgCMatrix"
#>         0a45497c6bfbfb22 2708240b908e2eba c3c393e9a17c1981 d4074c845bb62800
#> HLA-ABC              865             2077             2480             2994
#> B2M                 1182             3448             5307             9753
#> CD11b                929               12               17                6
#> CD11c                133             1354               26                8


## ------------------------------------------------
## Method `PixelDB$proximity`
## ------------------------------------------------

# Fetch the proximity scores
prox <- db$proximity()
prox %>% head()
#> # A tibble: 6 × 9
#>   marker_1 marker_2 join_count join_count_expected_mean join_count_expected_sd
#>   <chr>    <chr>         <dbl>                    <dbl>                  <dbl>
#> 1 CD56     CD56              0                     0                     0    
#> 2 CD56     mIgG2b            0                     0.03                  0.171
#> 3 CD56     CD71              0                     0.01                  0.1  
#> 4 CD56     CD6               0                     2.07                  1.58 
#> 5 CD56     Siglec-9          0                     0.03                  0.171
#> 6 CD56     CD79a             0                     0                     0    
#> # ℹ 4 more variables: join_count_z <dbl>, join_count_p <dbl>, component <chr>,
#> #   log2_ratio <dbl>


## ------------------------------------------------
## Method `PixelDB$cell_meta`
## ------------------------------------------------

# Fetch the cell meta data
db$cell_meta() %>% head()
#>                  n_umi1 n_umi2 n_edges n_antibodies reads_in_component n_umi
#> 0a45497c6bfbfb22  18560  24983   97014          149             289856 43543
#> 2708240b908e2eba  16301  21364   79638          158             229977 37665
#> c3c393e9a17c1981  21257  28094  110657          158             298516 49351
#> d4074c845bb62800  26891  38821  141153          156             372618 65712
#> efe0ed189cb499fc  19617  29019  100132          142             266795 48636
#>                  isotype_fraction intracellular_fraction tau_type       tau
#> 0a45497c6bfbfb22     0.0003674529                      0   normal 0.9865152
#> 2708240b908e2eba     0.0010885437                      0   normal 0.9563997
#> c3c393e9a17c1981     0.0022086685                      0   normal 0.9661265
#> d4074c845bb62800     0.0003956659                      0   normal 0.9818987
#> efe0ed189cb499fc     0.0002467308                      0   normal 0.9794434
#>                              sample antibodies average_k_core k_core_1 k_core_2
#> 0a45497c6bfbfb22 PNA055_Sample07_S7        149       2.569093    10260     6046
#> 2708240b908e2eba PNA055_Sample07_S7        158       2.466507     9596     6826
#> c3c393e9a17c1981 PNA055_Sample07_S7        158       2.639055    11845     8817
#> d4074c845bb62800 PNA055_Sample07_S7        156       2.474403    16575    11611
#> efe0ed189cb499fc PNA055_Sample07_S7        142       2.310058    12156     9552
#>                  k_core_3 k_core_4 k_core_5 svd_var_expl_s1 svd_var_expl_s2
#> 0a45497c6bfbfb22    19434     7803        0       0.4784136       0.1302096
#> 2708240b908e2eba    15319     5924        0       0.2922660       0.1960881
#> c3c393e9a17c1981    13995    14694        0       0.3671295       0.1443445
#> d4074c845bb62800    27303    10223        0       0.4026852       0.2052466
#> efe0ed189cb499fc    26620      308        0       0.3301883       0.1886183
#>                  svd_var_expl_s3 A_nodes_mean_degree B_nodes_mean_degree
#> 0a45497c6bfbfb22      0.08249094            5.227047            3.883201
#> 2708240b908e2eba      0.09845118            4.885467            3.727673
#> c3c393e9a17c1981      0.11940564            5.205673            3.938813
#> d4074c845bb62800      0.06080128            5.249080            3.635996
#> efe0ed189cb499fc      0.10618119            5.104348            3.450567


## ------------------------------------------------
## Method `PixelDB$protein_meta`
## ------------------------------------------------

# Fetch the protein meta data
db$protein_meta() %>% head()
#>         antibody_count antibody_pct components control nuclear
#> HLA-ABC         483065 0.0562623443        193      no      no
#> B2M            1095103 0.1275461109        193      no      no
#> CD11b             7907 0.0009209244        189      no      no
#> CD11c            34574 0.0040268169        192      no      no
#> CD18             90833 0.0105792751        192      no      no
#> CD82             87827 0.0102291677        193      no      no


## ------------------------------------------------
## Method `PixelDB$run_meta`
## ------------------------------------------------

# Fetch the run meta data
db$run_meta()
#> # A tibble: 1 × 6
#>   sample_name        version technology   panel_name panel_version post_analysis
#>   <chr>              <chr>   <chr>        <chr>      <chr>         <named list> 
#> 1 PNA055_Sample07_S7 0.1.0   single-cell… pna-rnd-1… 0.1.0         <named list> 


## ------------------------------------------------
## Method `PixelDB$components_edgelist`
## ------------------------------------------------

# Fetch edgelists
db$components_edgelist("0a45497c6bfbfb22") %>% head()
#> # A tibble: 6 × 5
#>   marker_1 marker_2    umi1    umi2 component       
#>   <chr>    <chr>    <int64> <int64> <chr>           
#> 1 HLA-ABC  HLA-ABC    6 e16   5 e16 0a45497c6bfbfb22
#> 2 HLA-ABC  HLA-ABC    6 e16   4 e16 0a45497c6bfbfb22
#> 3 HLA-ABC  HLA-ABC    1.e16   4 e15 0a45497c6bfbfb22
#> 4 HLA-ABC  HLA-ABC    5.e16   2.e16 0a45497c6bfbfb22
#> 5 HLA-ABC  HLA-ABC    6 e16   6 e16 0a45497c6bfbfb22
#> 6 HLA-ABC  HLA-ABC    1 e16   1.e16 0a45497c6bfbfb22


## ------------------------------------------------
## Method `PixelDB$components_layout`
## ------------------------------------------------

# Fetch layouts (NOTE: This will only work if the layouts exist in the database)
if (FALSE) { # \dontrun{
db$components_layout("0a45497c6bfbfb22")[[1]] %>% head()
} # }


## ------------------------------------------------
## Method `PixelDB$components_marker_counts`
## ------------------------------------------------

# Fetch marker counts
db$components_marker_counts("0a45497c6bfbfb22")[[1]][1:3, 1:4]
#> # A tibble: 3 × 4
#>   name                    CD16  CD44  CD59
#>   <chr>                  <dbl> <dbl> <dbl>
#> 1 55619162189985174-umi1     1     0     0
#> 2 71258604383722477-umi1     1     0     0
#> 3 5652424385946256-umi1      1     0     0


## ------------------------------------------------
## Method `PixelDB$export_parquet`
## ------------------------------------------------

# Export a table to a parquet file
tmp_parquet_file <- fs::file_temp(ext = "parquet")
db$export_parquet(tmp_parquet_file, "proximity")
fs::file_exists(tmp_parquet_file)
#> C:/Users/max/AppData/Local/Temp/Rtmpampkmn/file5bf44bcd64ee.parquet 
#>                                                                TRUE 


## ------------------------------------------------
## Method `PixelDB$close`
## ------------------------------------------------

# Close the connection when finished
db$close()
```
