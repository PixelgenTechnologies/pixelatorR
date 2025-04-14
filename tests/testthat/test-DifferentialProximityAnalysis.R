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
          n_tgt = c(30L, 30L),
          n_ref = c(10L, 10L),
          diff_median = c(1.12578767971233, -0.660561539890729),
          statistic = c(100, 200),
          auc = c(
            0.333333333333333,
            0.666666666666667
          ),
          p = c(0.11041579736434, 0.11041579736434),
          p_adj = c(1, 1),
          method = c(
            "Wilcoxon rank-sum test test",
            "Wilcoxon rank-sum test test"
          ),
          alternative = c(
            "two.sided",
            "two.sided"
          ),
          marker_1 = c("B2M", "B2M"),
          marker_2 = c(
            "B2M",
            "CD10"
          )
        ),
        row.names = c(NA, -2L),
        class = c(
          "tbl_df", "tbl",
          "data.frame"
        )
      )

    expect_equal(head(de_results1, 2), expected_data)

    # data.table backend
    expect_no_error({
      de_results2 <- DifferentialProximityAnalysis(
        seur_obj_big,
        contrast_column = "cell_type",
        targets = "CD4T",
        reference = "pDC",
        backend = "data.table",
        verbose = FALSE
      )
    })
    expect_equal(de_results1, de_results2)
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
