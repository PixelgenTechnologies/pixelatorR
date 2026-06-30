library(dplyr)

pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expect_no_error(seur_obj <- suppressWarnings(ReadPNA_Seurat(pxl_file, verbose = FALSE)))
  seur_obj$cell_type <- c("CD16+ Mono", "pDC", "CD4T", "CD4T", "CD4T")

  # Create a bigger Seurat object
  seur_obj_big <- merge(seur_obj, rep(list(seur_obj), 9), add.cell.ids = LETTERS[1:10])

  if (assay_version == "v5") {
    seur_obj_big <- seur_obj_big %>% JoinLayers()
  }

  test_that("DifferentialProximityAnalysis works as expected", {
    # Test DifferentialProximityAnalysis
    expect_no_error({
      de_results1 <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        diff_threshold = 1,
        verbose = FALSE
      )
    })

    # Test multiple targets
    expect_no_error({
      de_results_multiple_targets <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = c("CD4T", "CD16+ Mono"),
        reference = "pDC",
        diff_threshold = 1,
        verbose = FALSE
      )
    })
    expect_equal(
      structure(
        c(CD4T = 931L, `CD16+ Mono` = 1206L),
        dim = 2L,
        dimnames = list(
          . = c("CD4T", "CD16+ Mono")
        ),
        class = "table"
      ),
      de_results_multiple_targets$target %>% table() %>% sort()
    )

    # Test that contrast column can be factor
    prox_scores <-
      ProximityScores(seur_obj_big,
        meta_data_columns = "cell_type"
      ) %>%
      mutate(cell_type = factor(cell_type)) %>%
      filter(marker_1 == marker_2)

    expect_no_error({
      de_results <- DifferentialProximityAnalysis(
        prox_scores,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        verbose = FALSE
      )
    })

    expected_data <-
      structure(
        list(
          data_type = c("log2_ratio", "log2_ratio"),
          target = c("CD4T", "CD4T"),
          reference = c("pDC", "pDC"),
          pct_tgt = c(0, 0.333),
          pct_ref = c(0, 0),
          diff_median = c(1.13093086982645, 2.06004738366994),
          p = c(0.0924377451134652, 3.73436610215361e-07),
          p_adj = c(1, 0.000347669484110502),
          alternative = c("two.sided", "two.sided"),
          marker_1 = c("CD56", "CD56"),
          marker_2 = c("HLA-ABC", "HLA-DR-DP-DQ")
        ),
        row.names = c(NA, -2L),
        class = c("tbl_df", "tbl", "data.frame")
      )

    expect_equal(head(de_results1, 2), expected_data)

    # legacy method
    expect_no_error({
      de_results1 <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        method = "legacy",
        min_exp_join_count = 50,
        verbose = FALSE
      )
    })

    # data.table backend
    expect_no_error({
      de_results2 <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        method = "legacy",
        min_exp_join_count = 50,
        backend = "data.table",
        verbose = FALSE
      )
    })
    expect_equal(de_results1, de_results2)

    # min_cells_per_group option
    expect_no_error({
      de_results_filt <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        diff_threshold = 1,
        min_cells_per_group = 10,
        verbose = FALSE
      )
    })
    expect_equal(dim(de_results_filt), c(931, 11))
  })

  test_that("DifferentialProximityAnalysis fails with invalid input", {
    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          targets = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          group_vars = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          assay = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          proximity_metric = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          metric_type = "Invalid"
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          min_n_obs = -1
        )
      },
      class = "error"
    )

    expect_error(
      {
        DifferentialProximityAnalysis(
          seur_obj_big,
          contrast_column = "cell_type",
          reference = "B",
          p_adjust_method = "Invalid"
        )
      },
      class = "error"
    )
  })
}
