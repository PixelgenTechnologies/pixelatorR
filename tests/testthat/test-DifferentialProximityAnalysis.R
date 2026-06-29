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
        verbose = FALSE
      )
    })
    expect_equal(
      structure(
        c(CD4T = 773L, `CD16+ Mono` = 878L),
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
      mutate(cell_type = factor(cell_type))

    expect_no_error({
      DifferentialProximityAnalysis(
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
          data_type = c("join_count_z", "join_count_z"),
          target = c("CD4T", "CD4T"),
          reference = c("pDC", "pDC"),
          n_tgt = c(0, 0),
          n_ref = c(10, 10),
          diff_median = c(-0.325941157393986, -0.806975315413949),
          p = c(4.84076936963803e-10, 4.84076936963803e-10),
          p_adj = c(3.74191472273019e-07, 3.74191472273019e-07),
          alternative = c("two.sided", "two.sided"),
          pair = c("CD50:CD53", "CD29:CD53")
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
        backend = "data.table",
        verbose = FALSE
      )
    })
    expect_equal(de_results1, de_results2)

    # min_join_count option
    expect_no_error({
      de_results_filt <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        min_join_count = 50,
        verbose = FALSE
      )
    })
    expect_equal(dim(de_results_filt), c(268, 10))

    # min_cells_per_group option
    expect_no_error({
      de_results_filt <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        min_cells_per_group = 10,
        verbose = FALSE
      )
    })
    expect_equal(dim(de_results_filt), c(773, 10))
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
