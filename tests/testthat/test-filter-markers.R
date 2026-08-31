seur <- ReadPNA_Seurat(minimal_pna_pxl_file(), overwrite = TRUE, load_proximity_scores = FALSE, verbose = FALSE)
isotype_markers <- c("mIgG1", "mIgG2a", "mIgG2b")
seur$sample <- c("S1", "S1", "S2", "S2", "S2")

test_that("FilterMarkers works as expected", {
  expect_no_error(
    kept <- FilterMarkers(
      seur,
      isotype_markers = isotype_markers
    )
  )

  expect_equal(
    kept,
    c(
      "HLA-ABC", "B2M", "CD11b", "CD11c", "CD18", "CD82", "CD8",
      "TCRab", "HLA-DR", "CD45", "CD14", "CD16", "CD19", "CD45RB",
      "CD44", "CD52", "CD59", "CD45RA", "CD36", "CD2", "CD29", "CD5",
      "CD162", "CD27", "CD35", "CD26", "CD49D", "CD37", "CD41", "CD20",
      "CD64", "CD22", "CD274", "CD328", "CD25", "CD279", "CD335", "CD86",
      "CD13", "CD156c", "CD158b", "CD159a", "CD163", "CD169", "CD1a",
      "CD206", "CD226", "CD273", "CD28", "CD352", "CD45RO", "CD49e",
      "CD80", "CD81", "CD85j", "CD93", "IgD", "IgE", "IgM", "KLRG1",
      "CD10", "CD103", "CD199", "CD21", "CD277", "CD305", "CD371",
      "CD89", "CD90", "CD70", "HLA-DQ", "HLA-DR-DP-DQ", "Siglec-9",
      "CD95", "CD127", "CD141", "CD31", "CD6", "CD57", "CD73", "CD366",
      "CD357", "CD193", "CD319", "TIGIT", "CD102", "CD123", "CD150",
      "CD154", "CD161", "CD180", "CD191", "CD192", "CD1d", "CD200",
      "CD229", "CD24", "CD244", "CD268", "CD278", "CD32", "CD33", "CD38",
      "CD39", "CD40", "CD48", "CD50", "CD55", "CD58", "CD69", "CD72",
      "CD84", "CD9", "CD94", "TCRVB5", "CD3e", "CD4", "CD11a", "CD43",
      "CD7", "CD53", "CD302", "VISTA", "CD269", "CD66b", "CX3CR1",
      "CD369", "CD54", "CD71", "CD47", "CD117", "CD314"
    )
  )

  expect_false(any(isotype_markers %in% kept))

  expect_equal(
    FilterMarkers(
      seur,
      isotype_markers = isotype_markers,
      isotype_ratio = NULL,
      abundance_threshold = 2000
    ),
    c(
      "HLA-ABC", "B2M", "CD11b", "CD11c", "CD18", "CD82", "TCRab",
      "HLA-DR", "CD45", "CD16", "CD45RB", "CD44", "CD52", "CD59", "CD45RA",
      "CD36", "CD2", "CD29", "CD5", "CD162", "CD27", "CD35", "CD26",
      "CD49D", "CD37", "CD328", "CD13", "CD156c", "CD226", "CD273",
      "CD352", "CD45RO", "CD49e", "CD81", "CD85j", "CD93", "KLRG1",
      "CD199", "CD305", "CD371", "CD89", "HLA-DR-DP-DQ", "Siglec-9",
      "CD127", "CD31", "CD6", "CD319", "CD102", "CD123", "CD150", "CD180",
      "CD24", "CD244", "CD278", "CD32", "CD38", "CD39", "CD40", "CD48",
      "CD50", "CD55", "CD58", "CD84", "CD3e", "CD4", "CD11a", "CD43",
      "CD7", "CD53", "CD302", "VISTA", "CD66b", "CD54", "CD47"
    )
  )

  expect_equal(
    FilterMarkers(
      seur,
      isotype_markers = isotype_markers,
      isotype_ratio = 1.5,
      abundance_threshold = 50000
    ),
    c(
      "HLA-ABC", "B2M", "CD45", "CD16", "CD45RB", "CD44", "CD59",
      "HLA-DR-DP-DQ", "CD24", "CD3e", "CD43", "CD66b"
    )
  )

  stats <- FilterMarkers(
    seur,
    isotype_markers = isotype_markers,
    return_stats = TRUE
  )
  expect_equal(
    as.data.frame(dplyr::filter(stats, marker %in% c("HLA-ABC", "CD45", "mIgG1"))),
    data.frame(
      marker = c("HLA-ABC", "CD45", "mIgG1"),
      positive_fraction = c(1, 1, 0.4),
      isotype_median_cpm = c(252.623843097628, 252.623843097628, 252.623843097628),
      isotype_ratio = c(1.5, 1.5, 1.5),
      abundance_threshold = c(NA_real_, NA_real_, NA_real_),
      kept = c(TRUE, TRUE, FALSE),
      stringsAsFactors = FALSE
    )
  )

  expect_equal(
    FilterMarkers(
      seur,
      isotype_markers = isotype_markers,
      group_column = "sample"
    ),
    list(
      S1 = c(
        "HLA-ABC", "B2M", "CD11b", "CD11c", "CD18", "CD82",
        "CD8", "TCRab", "HLA-DR", "CD45", "CD16", "CD45RB", "CD44", "CD52",
        "CD59", "CD45RA", "CD36", "CD29", "CD162", "CD27", "CD35", "CD49D",
        "CD37", "CD41", "CD64", "CD328", "CD86", "CD13", "CD156c", "CD163",
        "CD169", "CD273", "CD352", "CD45RO", "CD49e", "CD81", "CD85j",
        "CD93", "IgM", "KLRG1", "CD10", "CD199", "CD21", "CD305", "CD371",
        "CD89", "HLA-DQ", "HLA-DR-DP-DQ", "Siglec-9", "CD95", "CD141",
        "CD31", "CD6", "CD366", "CD357", "CD193", "CD319", "TIGIT", "CD102",
        "CD123", "CD150", "CD180", "CD191", "CD192", "CD1d", "CD24",
        "CD244", "CD32", "CD33", "CD38", "CD39", "CD40", "CD48", "CD50",
        "CD55", "CD58", "CD69", "CD72", "CD84", "CD9", "CD3e", "CD4",
        "CD11a", "CD43", "CD7", "CD53", "CD302", "VISTA", "CD66b", "CX3CR1",
        "CD369", "CD54", "CD71", "CD47"
      ),
      S2 = c(
        "HLA-ABC", "B2M", "CD11b",
        "CD11c", "CD18", "CD82", "CD8", "TCRab", "HLA-DR", "CD45", "CD14",
        "CD16", "CD19", "CD45RB", "CD44", "CD52", "CD59", "CD45RA", "CD36",
        "CD2", "CD29", "CD5", "CD162", "CD27", "CD26", "CD49D", "CD37",
        "CD41", "CD20", "CD22", "CD274", "CD328", "CD25", "CD279", "CD335",
        "CD86", "CD13", "CD156c", "CD158a", "CD158b", "CD159a", "CD159c",
        "CD163", "CD169", "CD1a", "CD206", "CD226", "CD273", "CD28",
        "CD352", "CD45RO", "CD49e", "CD56", "CD80", "CD81", "CD85j",
        "CD93", "IgD", "IgE", "IgM", "KLRG1", "CD103", "CD199", "CD21",
        "CD277", "CD305", "CD371", "CD89", "CD90", "CD70", "HLA-DQ",
        "HLA-DR-DP-DQ", "Siglec-9", "TCRva7.2", "CD95", "TCRVg9", "CD127",
        "CD31", "CD6", "CD57", "CD73", "CD366", "CD357", "CD193", "CD319",
        "TIGIT", "CD102", "CD123", "CD150", "CD154", "CD161", "CD191",
        "CD192", "CD200", "CD229", "CD24", "CD244", "CD268", "CD278",
        "CD32", "CD33", "CD38", "CD40", "CD48", "CD50", "CD55", "CD58",
        "CD62P", "CD69", "CD72", "CD84", "CD9", "CD94", "TCRVB5", "CD3e",
        "CD4", "CD11a", "CD43", "CD7", "CD53", "CD302", "VISTA", "CD269",
        "CD66b", "CX3CR1", "CD209", "CD369", "CD54", "CD47", "CD117",
        "CD314"
      )
    )
  )

  expect_error(
    FilterMarkers(
      seur,
      isotype_markers = c("mIgG1", "misspelled")
    )
  )

  expect_error(
    FilterMarkers(
      seur[[]],
      isotype_markers = isotype_markers
    )
  )

  expect_error(
    FilterMarkers(
      seur,
      isotype_markers = isotype_markers,
      isotype_ratio = NULL,
      abundance_threshold = NULL
    )
  )

  expect_error(
    FilterMarkers(
      seur,
      isotype_markers = isotype_markers,
      group_column = "missing_column"
    )
  )
})
