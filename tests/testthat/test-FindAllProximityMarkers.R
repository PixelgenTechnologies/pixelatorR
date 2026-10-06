library(dplyr)

arrange_markers <- function(x) {
  x %>% arrange(marker_1, marker_2, target, reference)
}

matrix_one_vs_rest <- function(
  df,
  ident,
  group_by = "cluster",
  proximity_metric = "join_count_z",
  rest = "rest",
  min_cells_per_group = 4,
  diff_threshold = 0.01,
  min_pct = 0,
  min_diff_pct = -Inf,
  only_pos = FALSE
) {
  labels <- df %>% distinct(component, !!sym(group_by))
  label_chr <- as.character(labels[[group_by]])
  names(label_chr) <- as.character(labels$component)
  mat <- ProximityScoresToAssay(df, values_from = proximity_metric)
  missing <- setdiff(names(label_chr), colnames(mat))
  if (length(missing) > 0) {
    m_missing <- Matrix::rsparsematrix(
      nrow = nrow(mat),
      ncol = length(missing),
      density = 0
    )
    rownames(m_missing) <- rownames(mat)
    colnames(m_missing) <- missing
    mat <- cbind(mat, m_missing)
  }
  mat <- mat[, names(label_chr), drop = FALSE]
  cell_labels <- unname(label_chr[colnames(mat)])
  contrast <- ifelse(
    is.na(cell_labels),
    NA_character_,
    ifelse(cell_labels == ident, ident, rest)
  )
  group_data <- data.frame(
    contrast = contrast,
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )
  DifferentialProximityAnalysis(
    mat,
    group_data = group_data,
    contrast_column = "contrast",
    reference = rest,
    targets = ident,
    proximity_metric = proximity_metric,
    min_cells_per_group = min_cells_per_group,
    diff_threshold = diff_threshold,
    min_pct = min_pct,
    min_diff_pct = min_diff_pct,
    only_pos = only_pos,
    verbose = FALSE
  )
}

one_vs_rest_data <- function(n_per_group = 4L, groups = c("0", "1", "2")) {
  components <- paste0("cell", seq_len(n_per_group * length(groups)))
  clusters <- rep(groups, each = n_per_group)
  tidyr::expand_grid(
    marker_1 = c("A", "B"),
    marker_2 = c("A", "B"),
    component = components
  ) %>%
    mutate(
      cluster = clusters[match(component, components)],
      join_count_z = rnorm(n())
    )
}

test_that("FindAllProximityMarkers compares each level to the pooled rest", {
  set.seed(1)
  df <- one_vs_rest_data(n_per_group = 4L, groups = as.character(0:5))
  df <- df %>%
    mutate(
      join_count_z = if_else(
        marker_1 == "A" & marker_2 == "A" & cluster == "0",
        join_count_z + 5,
        join_count_z
      )
    )
  before <- df

  expect_no_error({
    res <- FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      min_cells_per_group = 4,
      verbose = FALSE
    )
  })
  expect_identical(df, before)
  expect_s3_class(res, "tbl_df")
  expect_equal(unique(res$target), as.character(0:5))
  expect_true(all(res$reference == "rest"))
  expect_false(any(res$target == res$reference))

  boosted <- res %>%
    filter(target == "0", marker_1 == "A", marker_2 == "A")
  expect_gt(boosted$diff_median, 0)

  for (ident in unique(res$target)) {
    ref <- matrix_one_vs_rest(df, ident = ident, min_cells_per_group = 4)
    got <- res %>% filter(target == ident)
    expect_equal(arrange_markers(got), arrange_markers(ref))
  }
})

test_that("FindAllProximityMarkers respects factor order, idents, and the rest label", {
  set.seed(1)
  df <- one_vs_rest_data()
  df$cluster <- factor(df$cluster, levels = c("2", "0", "1"))

  res <- FindAllProximityMarkers(
    df,
    group_by = "cluster",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    verbose = FALSE
  )
  expect_equal(unique(res$target), c("2", "0", "1"))

  res_one <- FindAllProximityMarkers(
    df,
    group_by = "cluster",
    idents = "0",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    verbose = FALSE
  )
  expect_equal(unique(res_one$target), "0")

  res_numeric <- FindAllProximityMarkers(
    df %>% mutate(cluster = as.character(cluster)),
    group_by = "cluster",
    idents = 1,
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    verbose = FALSE
  )
  expect_equal(unique(res_numeric$target), "1")

  colliding <- df %>% mutate(cluster = if_else(as.character(cluster) == "2", "rest", as.character(cluster)))
  res_rest <- FindAllProximityMarkers(
    colliding,
    group_by = "cluster",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    verbose = FALSE
  )
  expect_true(all(res_rest$reference == "rest_"))
  expect_setequal(unique(res_rest$target), c("0", "1", "rest"))
})

test_that("FindAllProximityMarkers skips groups that are too small", {
  set.seed(1)
  df <- bind_rows(
    one_vs_rest_data(n_per_group = 4L, groups = c("0", "1")),
    one_vs_rest_data(n_per_group = 2L, groups = "2") %>%
      mutate(component = paste0("small", component))
  )

  expect_warning(
    res <- FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      min_cells_per_group = 4,
      verbose = FALSE
    ),
    "fewer than"
  )
  expect_setequal(unique(res$target), c("0", "1"))

  expect_warning(
    expect_error(
      FindAllProximityMarkers(
        df,
        group_by = "cluster",
        idents = "2",
        proximity_metric = "join_count_z",
        min_cells_per_group = 4,
        verbose = FALSE
      ),
      "No one-versus-rest"
    ),
    "fewer than"
  )
})

test_that("matrix thresholds are applied once to the shared proximity matrix", {
  set.seed(1)
  df <- one_vs_rest_data(n_per_group = 4L, groups = as.character(0:5))
  df <- df %>%
    mutate(
      join_count_z = if_else(
        marker_1 == "A" & marker_2 == "A" & cluster == "0",
        join_count_z + 5,
        join_count_z
      )
    )

  loose <- FindAllProximityMarkers(
    df,
    group_by = "cluster",
    idents = "0",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    diff_threshold = 0,
    min_pct = 0,
    min_diff_pct = -Inf,
    verbose = FALSE
  )
  strict <- FindAllProximityMarkers(
    df,
    group_by = "cluster",
    idents = "0",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    diff_threshold = 1,
    min_pct = 0,
    min_diff_pct = -Inf,
    verbose = FALSE
  )
  expect_lt(nrow(strict), nrow(loose))
  expect_true(all(abs(strict$diff_median) >= 1))
  ref <- matrix_one_vs_rest(
    df,
    ident = "0",
    min_cells_per_group = 4,
    diff_threshold = 1
  )
  expect_equal(arrange_markers(strict), arrange_markers(ref))
})

test_that("only_pos tests pairs with a higher median in the target", {
  set.seed(1)
  df <- one_vs_rest_data(n_per_group = 4L, groups = c("0", "1", "2"))
  df <- df %>%
    mutate(
      join_count_z = if_else(
        marker_1 == "A" & marker_2 == "A" & cluster == "0",
        join_count_z + 5,
        join_count_z
      )
    )

  both <- FindAllProximityMarkers(
    df,
    group_by = "cluster",
    idents = "0",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    diff_threshold = 0,
    only_pos = FALSE,
    verbose = FALSE
  )
  pos <- FindAllProximityMarkers(
    df,
    group_by = "cluster",
    idents = "0",
    proximity_metric = "join_count_z",
    min_cells_per_group = 4,
    diff_threshold = 0,
    only_pos = TRUE,
    verbose = FALSE
  )

  expect_true(any(both$diff_median < 0))
  expect_true(all(pos$diff_median > 0))
  expect_lt(nrow(pos), nrow(both))
  positive_pairs <- both %>%
    filter(diff_median > 0) %>%
    mutate(pair = paste(marker_1, marker_2))
  expect_setequal(paste(pos$marker_1, pos$marker_2), positive_pairs$pair)

  ref <- matrix_one_vs_rest(
    df,
    ident = "0",
    min_cells_per_group = 4,
    diff_threshold = 0,
    only_pos = TRUE
  )
  expect_equal(arrange_markers(pos), arrange_markers(ref))
})

test_that("FindAllProximityMarkers fails with invalid input", {
  set.seed(1)
  df <- one_vs_rest_data()

  expect_error(
    FindAllProximityMarkers(df),
    "group_by"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = c("cluster", "component"),
      proximity_metric = "join_count_z"
    ),
    "single string"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "missing",
      proximity_metric = "join_count_z"
    ),
    "missing"
  )
  expect_error(
    FindAllProximityMarkers(
      df %>% mutate(cluster = 1:n()),
      group_by = "cluster",
      proximity_metric = "join_count_z"
    ),
    "must be one of"
  )
  expect_error(
    FindAllProximityMarkers(
      df %>% mutate(cluster = "only"),
      group_by = "cluster",
      proximity_metric = "join_count_z",
      min_cells_per_group = 1
    ),
    "at least 2 groups"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      idents = "missing",
      proximity_metric = "join_count_z",
      min_cells_per_group = 4
    ),
    "Not all"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      idents = c("0", "0"),
      proximity_metric = "join_count_z",
      min_cells_per_group = 4
    ),
    "unique"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      reference = "rest",
      min_cells_per_group = 4
    ),
    "Cannot pass"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      group_vars = "marker_1",
      min_cells_per_group = 4
    ),
    "Cannot pass"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "marker_1",
      proximity_metric = "join_count_z",
      min_cells_per_group = 4
    ),
    "grouping column"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      min_cells_per_group = -1
    ),
    "non-negative"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      method = "legacy",
      min_cells_per_group = 4
    ),
    "method"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      min_pct = 2,
      min_cells_per_group = 4
    ),
    "min_pct"
  )
  expect_error(
    FindAllProximityMarkers(
      df,
      group_by = "cluster",
      proximity_metric = "join_count_z",
      only_pos = "yes",
      min_cells_per_group = 4
    ),
    "only_pos"
  )
})

pxl_file <- minimal_pna_pxl_file()

for (assay_version in c("v3", "v5")) {
  options(Seurat.object.assay.version = assay_version)

  expect_no_error(seur_obj <- suppressWarnings(ReadPNA_Seurat(pxl_file, verbose = FALSE)))
  seur_obj$cell_type <- c("CD16+ Mono", "pDC", "CD4T", "CD4T", "CD4T")
  seur_obj_big <- merge(seur_obj, rep(list(seur_obj), 9), add.cell.ids = LETTERS[1:10])

  if (assay_version == "v5") {
    seur_obj_big <- seur_obj_big %>% JoinLayers()
  }

  test_that(paste0("FindAllProximityMarkers works on a Seurat object (", assay_version, ")"), {
    meta_before <- colnames(seur_obj_big[[]])
    expect_no_error({
      res <- FindAllProximityMarkers(
        seur_obj_big,
        group_by = "cell_type",
        idents = "CD4T",
        metric_type = "self",
        diff_threshold = 0,
        min_cells_per_group = 10,
        verbose = FALSE
      )
    })
    expect_equal(colnames(seur_obj_big[[]]), meta_before)
    expect_true(all(res$target == "CD4T"))
    expect_true(all(res$reference == "rest"))
    expect_gt(nrow(res), 0)

    se_manual <- seur_obj_big
    se_manual$pxl_contrast <- ifelse(
      as.character(se_manual$cell_type) == "CD4T",
      "CD4T",
      "rest"
    )
    ref <- DifferentialProximityAnalysis(
      se_manual,
      contrast_column = "pxl_contrast",
      reference = "rest",
      targets = "CD4T",
      metric_type = "self",
      diff_threshold = 0,
      min_cells_per_group = 10,
      verbose = FALSE
    )
    expect_equal(arrange_markers(res), arrange_markers(ref))

    expect_no_error({
      res_all <- FindAllProximityMarkers(
        seur_obj_big,
        group_by = "cell_type",
        metric_type = "self",
        diff_threshold = 1,
        min_pct = 0,
        min_diff_pct = -Inf,
        min_cells_per_group = 10,
        verbose = FALSE
      )
    })
    expect_setequal(unique(res_all$target), c("CD16+ Mono", "pDC", "CD4T"))

    res_small <- suppressWarnings(FindAllProximityMarkers(
      seur_obj_big,
      group_by = "cell_type",
      metric_type = "self",
      diff_threshold = 1,
      min_cells_per_group = 11,
      verbose = FALSE
    ))
    expect_equal(unique(res_small$target), "CD4T")
  })

  test_that(paste0("FindAllProximityMarkers rejects invalid Seurat input (", assay_version, ")"), {
    expect_error(
      FindAllProximityMarkers(
        seur_obj_big,
        group_by = "missing"
      ),
      "missing"
    )
    expect_error(
      FindAllProximityMarkers(
        seur_obj_big,
        group_by = "cell_type",
        idents = "missing",
        verbose = FALSE
      ),
      "Not all"
    )
    expect_warning(
      expect_error(
        FindAllProximityMarkers(
          seur_obj_big,
          group_by = "cell_type",
          idents = "pDC",
          min_cells_per_group = 11,
          metric_type = "self",
          verbose = FALSE
        ),
        "No one-versus-rest"
      ),
      "fewer than"
    )
    expect_error(
      FindAllProximityMarkers(
        seur_obj_big,
        group_by = "cell_type",
        targets = "CD4T",
        verbose = FALSE
      ),
      "Cannot pass"
    )
  })

  test_that(paste0("FindAllProximityMarkers fetches unloaded proximity scores when lazy = TRUE (", assay_version, ")"), {
    seur_missing <- suppressWarnings(ReadPNA_Seurat(
      pxl_file,
      load_proximity_scores = FALSE,
      verbose = FALSE
    ))
    seur_missing$cell_type <- seur_obj$cell_type

    expect_error(
      FindAllProximityMarkers(
        seur_missing,
        group_by = "cell_type",
        metric_type = "self",
        min_cells_per_group = 2,
        verbose = FALSE
      ),
      "lazy = TRUE"
    )

    res_lazy <- FindAllProximityMarkers(
      seur_missing,
      group_by = "cell_type",
      idents = "CD4T",
      lazy = TRUE,
      metric_type = "self",
      diff_threshold = 0,
      min_cells_per_group = 2,
      verbose = FALSE
    )
    res_loaded <- FindAllProximityMarkers(
      seur_obj,
      group_by = "cell_type",
      idents = "CD4T",
      lazy = FALSE,
      metric_type = "self",
      diff_threshold = 0,
      min_cells_per_group = 2,
      verbose = FALSE
    )
    expect_equal(arrange_markers(res_lazy), arrange_markers(res_loaded))
  })
}
