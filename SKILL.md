# Skill: Generate Unit Tests for pixelatorR

## Objective

Generate `testthat` unit test files for R functions in the `pixelatorR` package.
Each exported or internal function should have a dedicated test script in `tests/testthat/` following the conventions described below. Closely related functions/methods (for example a function that has methods for different input classes) may be grouped together in the same test script if appropriate.

---

## Project Context

- **Package**: `pixelatorR` — an R/Seurat extension for Molecular Pixelation (MPX) and Proximity Network Assay (PNA) single-cell spatial proteomics data.
- **Test framework**: [`testthat`](https://testthat.r-lib.org/) (>= 3.0.0).
- **Test runner entry point**: `tests/testthat.R` (loads `testthat` and `pixelatorR`).
- **Source code**: `R/` directory. Each file typically contains one or a few closely related exported functions.

---

## Test File Conventions

### Naming

Each function gets **one** test file:

```
tests/testthat/test-<function_name>.R
```

Use the function name exactly (case-sensitive). Examples:

| Function              | Test file                              |
|-----------------------|----------------------------------------|
| `ColocalizationHeatmap` | `test-ColocalizationHeatmap.R`       |
| `ReadPNA_proximity`     | `test-ReadPNA_proximity.R`           |
| `local_G`               | `test-local_G.R`                     |
| `CalculateDispersion`   | `test-CalcDispersion.R`              |

### File Structure

Every test file follows this three-part layout:

```r
# ── 1. Setup: Load data ──────────────────────────────────
library(dplyr)
pxl_file <- minimal_pna_pxl_file()
prox <- ReadPNA_proximity(pxl_file)

# ── 2. Tests: Expected behavior ──────────────────────────
test_that("<function_name> works as expected", {
  expect_no_error(result <- <function_name>(...))
  # Prefer data-level assertions (see "Assertion Preferences" below)
  expect_equal(result[1:5, 1:3], <expected_snippet>)
})

# ── 3. Tests: Invalid input ──────────────────────────────
test_that("<function_name> fails with invalid input", {
  expect_error(<function_name>("Invalid"))
  expect_error(<function_name>(data, bad_param = "Invalid"))
})
```

### Key Rules

1. **One function per file.** Do not combine tests for unrelated functions.
2. **Top-level setup only.** Data loading and object creation go *outside* `test_that()` blocks so they are shared by all tests in the file.
3. **No `library(pixelatorR)`.** The package is already loaded by the test runner (`tests/testthat.R`).
4. **Only load extra libraries when needed.** Common additions: `library(dplyr)`, `library(Seurat)`. Only add them if the test code uses functions from those packages.

---

## Test Data

The package ships two minimal PXL files. Use the helpers below to get paths — they are always available after `library(pixelatorR)`. MPX will be deprecated in the future, so only PNA test data should be used for new tests.

| Helper                    | Data type | Cells | Description                    |
|---------------------------|-----------|-------|--------------------------------|
| `minimal_mpx_pxl_file()`  | MPX       | 5     | 5 immune cells (MPX protocol)  |
| `minimal_pna_pxl_file()`  | PNA       | 5     | 5 PBMC cells (PNA protocol)    |


### Loading Test Data

Choose the reader that fits the function under test:

```r
# Seurat object (PNA)
seur <- ReadPNA_Seurat(minimal_pna_pxl_file(), overwrite = TRUE,
                       load_proximity_scores = FALSE, verbose = FALSE)

# Proximity scores table
prox <- ReadPNA_proximity(minimal_pna_pxl_file())

# Count matrix
counts <- ReadPNA_counts(minimal_pna_pxl_file())

# Edgelist
edgelist <- ReadPNA_edgelist(minimal_pna_pxl_file())
```

When testing Seurat methods that require merged data:

```r
seur1 <- seur2 <- ReadPNA_Seurat(pxl_file, overwrite = TRUE)
seur1$sample <- "Sample1"
seur2$sample <- "Sample2"
seur_merged <- merge(seur1, seur2, add.cell.ids = c("Sample1", "Sample2"))
```

When testing CellGraph-level functions:

```r
seur <- ReadPNA_Seurat(minimal_pna_pxl_file())
seur <- LoadCellGraphs(seur, cells = colnames(seur)[1])
cg <- CellGraphs(seur)[[1]]
g <- CellGraphData(cg, slot = "cellgraph")
counts <- CellGraphData(cg, slot = "counts")
```

### Dealing with the `Seurat.object.assay.version` Option

The test runner sets `options(Seurat.object.assay.version = "v3")` globally.
If a test needs v5 behavior, set the option locally and restore it:

```r
test_that("function works with Assay5", {
  options(Seurat.object.assay.version = "v5")
  on.exit(options(Seurat.object.assay.version = "v3"), add = TRUE)
  # ... test code ...
})
```

---

## Assertion Preferences

Tests should compare actual output data — not just check dimensions or classes.
Rank from **most preferred** to **least preferred**:

### 1. Small Data Snippets (Preferred)

Compare a small slice of output to a hardcoded expected value. This catches regressions in the actual computation.

```r
# Compare first few rows/columns of a matrix result
expect_equal(result[1:5, 1:3], expected_snippet)

# Compare named vector values
expect_equal(
  result$scores,
  c(cell_a = -6.698, cell_b = 2.590, cell_c = 13.230),
  tolerance = 1e-3
)

# Compare a few rows of a data.frame
expect_equal(result[1:2, ], expected_df, tolerance = 1e-4)
```

### 2. Structural + Partial Data

When full data comparison is impractical, combine structural checks with partial data:

```r
expect_s3_class(result, "tbl_df")
expect_equal(names(result), c("marker_1", "marker_2", "estimate", "p_adj"))
expect_equal(result$marker_1[1:3], c("CD3E", "CD4", "CD8"))
```

### 3. Pure Structural (Use Sparingly)

Only as a supplement, never as the sole test:

```r
expect_identical(dim(result), c(100L, 5L))
expect_type(result, "double")
```

### How to Obtain Expected Values

Run the function in an interactive R session with the test data and capture a small representative output:

```r
result <- my_function(test_data)
dput(result[1:5, 1:3])  # Copy this into the test as the expected value
```

Use `dput()` to get exact R representations. For numeric values, keep enough decimal places to be meaningful but use `tolerance` in `expect_equal()` where floating-point differences are expected.

---

## Testing Invalid Input

For each parameter of the function under test, provide at least one invalid value that should trigger an error.

### Common Patterns

```r
test_that("<function_name> fails with invalid input", {
  # Wrong type for first argument
  expect_error(<function_name>("Invalid"))

  # Wrong column name
  expect_error(<function_name>(data, col_param = "NonexistentColumn"))

  # Wrong type for a boolean parameter
  expect_error(<function_name>(data, flag_param = "Invalid"))

  # Wrong type for a numeric parameter
  expect_error(<function_name>(data, num_param = "Invalid"))

  # Invalid method/choice

  expect_error(<function_name>(data, method = "Invalid"))

  # Missing required argument (if applicable)
  expect_error(<function_name>(data),
    'argument "required_param" is missing, with no default'
  )
})
```

### Tips

- Look at the function's parameter validation code (typically `assert_*` calls or `match.arg` at the top of the function body) to identify which values will be rejected.
- Use the string `"Invalid"` as the canonical bad input for character parameters — this is a project convention.
- For parameters with `match.arg`, passing an unrecognized string triggers the error.
- For boolean parameters, pass `"Invalid"` (a string instead of logical).
- For numeric parameters, pass `"Invalid"` (a string instead of numeric).

---

## Handling Warnings and Verbose Output

Some functions produce expected warnings or verbose messages. Suppress these to keep test output clean:

```r
expect_no_error(suppressWarnings(result <- noisy_function(data)))
```

---

## S4 Method Testing

Some functions are S4 generics with methods for different classes (`data.frame`, `Seurat`, `CellGraphAssay`, etc.).
Test each dispatched method in a separate `test_that()` block:

```r
test_that("RunDPA works as expected on a data.frame", { ... })
test_that("RunDPA works as expected on a Seurat object", { ... })
```

---

## Edge Cases Worth Testing

Consider adding tests for:

- **Zero/empty input**: empty data frames, zero-row slices
- **Single-element input**: one cell, one marker
- **Boundary conditions**: k=1 for neighborhood functions, resolution=0 for clustering
- **Factor vs character input**: some functions behave differently
- **Merged vs single-sample Seurat objects**
- **Division by zero**: functions operating on counts (e.g., normalization with all-zero rows)

---

## Parallel Processing Tests

If the function supports parallel execution (e.g., `cl` parameter), add a separate block:

```r
if (TRUE) skip("Skipping parallel processing tests")

test_that("<function_name> can be parallelized", {
  result_seq <- <function_name>(data, cl = 1)
  result_par <- <function_name>(data, cl = 2)
  expect_equal(result_seq, result_par)
})
```

Note: parallel tests are typically skipped by default (the `if (TRUE) skip(...)` pattern) since they can be flaky in CI.

---

## Step-by-Step Workflow

When asked to generate a test file for function `foo`:

1. **Read the source code** in `R/` to understand:
   - What the function does
   - Its parameters and their types/defaults
   - What validation/assertions it performs on inputs
   - What it returns (class, structure)
   - Whether it has S4 methods for multiple classes

2. **Identify the correct test data** (for example, Seurat vs raw table).

3. **Run the function** interactively to capture a small expected output snippet using `dput()`.

4. **Write the test file** following the three-part structure:
   - Setup (data loading, outside `test_that`)
   - Expected-behavior tests (with data-level comparisons)
   - Invalid-input tests (one per parameter where feasible)

5. **Run the tests** to confirm they pass:
   ```bash
   Rscript -e 'testthat::test_file("tests/testthat/test-foo.R")'
   ```

---

## Complete Example

Here is a full example for `ColocalizationHeatmap`:

```r
library(dplyr)
prox <- ReadPNA_proximity(minimal_pna_pxl_file())
prox_summarized <- prox %>%
  slice_sample(n = 1e4) %>%
  group_by(marker_1, marker_2) %>%
  summarize(mean_log2_ratio = mean(log2_ratio), .groups = "drop") %>%
  mutate(test1 = marker_1, test2 = marker_2, estimate = mean_log2_ratio) %>%
  ungroup()

test_that("ColocalizationHeatmap works as expected", {
  # Default method
  expect_no_error(ColocalizationHeatmap(prox_summarized))

  # Dots method
  expect_no_error(p <- ColocalizationHeatmap(prox_summarized, type = "dots",
                                              size_col = "mean_log2_ratio"))
  expect_true(is.factor(p$data$marker_1))
  expect_true(is.factor(p$data$marker_2))
  expect_equal(p$data$marker_1 %>% levels(), p$data$marker_2 %>% levels())

  # Custom marker columns
  expect_no_error(ColocalizationHeatmap(prox_summarized,
                                         marker1_col = "test1",
                                         marker2_col = "test2"))

  # Custom value column
  expect_no_error(ColocalizationHeatmap(prox_summarized,
                                         value_col = "mean_log2_ratio"))

  # Symmetrise FALSE
  expect_no_error(p <- ColocalizationHeatmap(prox_summarized, symmetrise = FALSE,
                                              type = "dots",
                                              size_col = "mean_log2_ratio"))
  expect_true(is.factor(p$data$marker_1))
  expect_true(is.factor(p$data$marker_2))
  expect_true(!identical(p$data$marker_1 %>% levels(),
                         p$data$marker_2 %>% levels()))
})

test_that("ColocalizationHeatmap fails with invalid input", {
  expect_error(ColocalizationHeatmap("Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, marker1_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, marker2_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, value_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, type = "dots",
                                      size_col = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized,
                                      size_col_transform = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, size_range = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, colors = FALSE))
  expect_error(ColocalizationHeatmap(prox_summarized, cluster_rows = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, cluster_cols = "Invalid"))
  expect_error(ColocalizationHeatmap(prox_summarized, type = "Invalid"))
})
```

---

## Checklist

Before finalizing a test file, verify:

- [ ] File is named `test-<function_name>.R`
- [ ] Data setup is outside `test_that()` blocks
- [ ] At least one `test_that(... "works as expected" ...)` block exists
- [ ] At least one `test_that(... "fails with invalid input" ...)` block exists
- [ ] Tests compare actual output data (not just dimensions/classes)
- [ ] Expected values were obtained by running the function with test data
- [ ] Each invalid-input test uses `expect_error()`
- [ ] No `library(pixelatorR)` call in the test file
- [ ] Extra `library()` calls only for packages actually used (e.g., `dplyr`, `Seurat`)
- [ ] Tests pass when run with `testthat::test_file()`
