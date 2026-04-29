# How to Use the Test Generation Skill

This guide explains how to use `SKILL.md` to automate the generation of `testthat` test files for the `pixelatorR` package using an LLM assistant (Gemini, Claude, ChatGPT, etc.).

---

## Quick Start

1. Open your preferred LLM assistant (Gemini, Claude, ChatGPT, Cursor, etc.)
2. Provide `SKILL.md` as context (paste it, attach it, or reference it)
3. Ask the LLM to generate a test file for a specific function

### Example Prompt

> Using the conventions in SKILL.md, generate a test file for the `FilterProximityScores` function. Here is the source code:
>
> *(paste the contents of `R/filter_and_summarize_proximity_scores.R`)*

---

## What to Provide as Context

For best results, give the LLM:

| Context                     | Why                                                       | Required? |
|-----------------------------|-----------------------------------------------------------|-----------|
| `SKILL.md`                  | Defines all conventions, patterns, and rules              | **Yes**   |
| Source code of the function | So the LLM understands parameters, validation, and return types | **Yes** |
| An existing test file       | As a concrete example of the style                        | Optional  |
| Roxygen documentation       | Extra detail on expected behavior                         | Optional  |

### Minimal Prompt Template

```
Using the test generation conventions from the attached SKILL.md,
write a test file for the function `<FUNCTION_NAME>`.

Source code:
<paste R/<source_file>.R>
```

### Richer Prompt Template

```
Using the test generation conventions from the attached SKILL.md,
write a test file for the function `<FUNCTION_NAME>`.

Source code:
<paste R/<source_file>.R>

Here is an existing test file for reference:
<paste tests/testthat/test-<similar_function>.R>

Please:
1. Use the PNA test data (minimal_pna_pxl_file())
2. Include tests for all S4 methods
3. Capture small data snippets as expected values using dput()
```

---

## After Generation

Once the LLM produces a test file:

### 1. Save the file

Save it to:
```
tests/testthat/test-<function_name>.R
```

### 2. Generate expected values

The LLM will likely use placeholder values for `expect_equal()` comparisons (since it cannot run R code). Replace them with real values:

```r
# In an interactive R session:
library(pixelatorR)
library(dplyr)

# Run the setup code from the generated test file
# then capture the output:
result <- my_function(test_data)
dput(result[1:5, 1:3])
# Copy the dput() output into the test file as the expected value
```

> **Tip:** If you're using an agent that _can_ execute R code (e.g., Cursor with a terminal), ask it to run the function and capture `dput()` output directly. This eliminates the manual step.

### 3. Run the tests

```bash
# Single file
Rscript -e 'testthat::test_file("tests/testthat/test-my_function.R")'

# All tests
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### 4. Fix any failures

Common issues:
- **Wrong expected values**: re-run with `dput()` and update
- **Missing library**: add `library(dplyr)` or similar at the top
- **Seurat version mismatch**: ensure `options(Seurat.object.assay.version = "v3")` is set