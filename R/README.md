# Interaction priors + ColocalizationHeatmap (prototype)

Package-shaped R code intended to replace `pixelatorR::ColocalizationHeatmap`
with optional cell-border highlighting. Database access is separate from plotting.

## Files

| File | Role |
|------|------|
| `ColocalizationHeatmap.R` | Heatmap (`dots` / `tiles`) with `highlight_pairs` cell borders |
| `ExtractPanelInteractions.R` | Panel × DB interaction extraction + filters |
| `interaction_priors.R` | Load/save slim UniProt edge priors |
| `BuildInteractionPriors.R` | `BuildStringPrior()`, `BuildBiogridPrior()`, … |

Licensing notes: [`docs/interaction_priors_LICENSING.md`](../docs/interaction_priors_LICENSING.md).

## Quick start

```r
# from project root
for (f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(f)

# once: build priors from results/db_cache/* releases
# source("scripts/build_all_priors.R")

highlight_pairs <- ExtractPanelInteractions(
  markers = unique(c(coloc_summary$marker_1, coloc_summary$marker_2)),
  database = "string",
  marker_uniprot_map = marker_map,
  score_threshold = 400,
  string_network = "physical"
) |>
  dplyr::select(marker_1, marker_2)

ColocalizationHeatmap(
  coloc_summary,
  value_col = "mean_log2_ratio",
  size_col = "pct_detected",
  type = "dots",
  highlight_pairs = highlight_pairs
)
```

Demo: `Rscript scripts/demo_coloc_db_highlight.R`
