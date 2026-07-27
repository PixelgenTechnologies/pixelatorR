# Canonical edge schema: uniprot_a, uniprot_b (undirected, a <= b), optional
# native score columns, evidence, resource, resource_version.

#' Default cache directory for raw interaction-database dumps
#'
#' Under \code{tools::R_user_dir("pixelatorR", "cache")/db_raw}. Builders
#' download missing source files here before writing slim RDS caches.
#'
#' @return Character path.
#' @noRd
.default_raw_cache <- function() {
  file.path(tools::R_user_dir("pixelatorR", which = "cache"), "db_raw")
}

#' Path to a slim interaction-database RDS file
#'
#' @param database Database key (e.g. \code{"string"}).
#' @param version Version label used in the filename (e.g. \code{"latest"}).
#' @param cache_dir Directory containing slim RDS caches.
#'
#' @return Character path \code{cache_dir/<database>_<version>.rds}.
#' @noRd
.interaction_database_rds_path <- function(database, version, cache_dir) {
  file.path(cache_dir, paste0(database, "_", version, ".rds"))
}

#' Download a file if it is missing or empty
#'
#' Raises \code{options(timeout)} for the transfer (see \code{?download.file})
#' and restores the previous value on exit. Uses \code{mode = "wb"} and shows
#' a progress bar (\code{quiet = FALSE}).
#'
#' @param url Remote URL to download.
#' @param dest_file Local destination path.
#' @param min_timeout Minimum timeout in seconds. Defaults to 300; large dumps
#'   (e.g. AlphaFold heterodimer metadata) should pass a higher floor.
#'
#' @return \code{dest_file} (existing or freshly downloaded).
#' @noRd
.download_if_missing <- function(url, dest_file, min_timeout = 300) {
  if (file.exists(dest_file) && isTRUE(file.info(dest_file)$size > 0)) {
    return(dest_file)
  }
  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
  old_timeout <- getOption("timeout")
  # Write to a sibling temp file first so a failed download cannot leave a
  # truncated dest_file that later runs would treat as a valid cache.
  tmp_file <- tempfile(
    pattern = paste0(".", basename(dest_file), "."),
    tmpdir = dirname(dest_file),
    fileext = ".part"
  )
  on.exit(
    {
      options(timeout = old_timeout)
      if (file.exists(tmp_file)) {
        unlink(tmp_file)
      }
    },
    add = FALSE
  )
  options(timeout = max(min_timeout, old_timeout))
  cli::cli_inform(c("i" = "Downloading {.url {url}}"))
  status <- tryCatch(
    utils::download.file(url, destfile = tmp_file, mode = "wb", quiet = FALSE),
    error = function(e) {
      cli::cli_abort(c(
        "x" = "Failed to download {.url {url}}",
        "i" = conditionMessage(e)
      ))
    }
  )
  if (!identical(as.integer(status), 0L) ||
    !file.exists(tmp_file) ||
    !isTRUE(file.info(tmp_file)$size > 0)) {
    cli::cli_abort("Download failed or produced an empty file: {.url {url}}")
  }
  if (!isTRUE(file.rename(tmp_file, dest_file))) {
    if (!isTRUE(file.copy(tmp_file, dest_file, overwrite = TRUE))) {
      cli::cli_abort("Downloaded {.url {url}} but failed to move it to {.path {dest_file}}")
    }
    unlink(tmp_file)
  }
  return(dest_file)
}

#' First matching column name among candidates (exact or regex)
#'
#' @param df Data frame or tibble.
#' @param candidates Character vector of exact names, or regex patterns when
#'   \code{regex = TRUE}.
#' @param regex If \code{TRUE}, treat each candidate as a case-insensitive
#'   regex against \code{names(df)}.
#'
#' @return A single column name, or \code{NA_character_} if none match.
#' @noRd
.first_present_col <- function(df, candidates, regex = FALSE) {
  nms <- names(df)
  if (isTRUE(regex)) {
    for (pat in candidates) {
      hits <- grep(pat, nms, value = TRUE, ignore.case = TRUE)
      if (length(hits) > 0) {
        return(hits[[1]])
      }
    }
    return(NA_character_)
  }
  hit <- intersect(candidates, nms)
  if (length(hit) < 1) {
    return(NA_character_)
  }
  return(hit[[1]])
}

#' Validate a named numeric score-bound vector
#'
#' @param x Named numeric vector or NULL.
#' @param arg_name Argument name for errors.
#' @param available Available score column names.
#'
#' @return \code{x} unchanged, or \code{NULL}.
#' @noRd
.validate_named_score_bounds <- function(x, arg_name, available) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) < 1) {
    cli::cli_abort("{.arg {arg_name}} must be a fully named numeric vector.")
  }
  nms <- names(x)
  if (is.null(nms) || anyNA(nms) || any(nms == "")) {
    cli::cli_abort("{.arg {arg_name}} must be a fully named numeric vector.")
  }
  if (anyDuplicated(nms)) {
    cli::cli_abort("{.arg {arg_name}} must not contain duplicate names.")
  }
  if (length(available) < 1) {
    cli::cli_abort(c(
      "x" = "This database has no score columns.",
      "i" = "Cannot use {.arg {arg_name}}."
    ))
  }
  unknown <- setdiff(nms, available)
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "x" = "Unknown score column{?s} in {.arg {arg_name}}: {.val {unknown}}.",
      "i" = "Available score columns: {.val {available}}."
    ))
  }
  return(x)
}

#' Filter interaction edges by named score bounds
#'
#' @param edges Edge tibble.
#' @param score_min Named numeric minima (\code{>=}), or NULL.
#' @param score_max Named numeric maxima (\code{<}), or NULL.
#' @param score_combine \code{"any"} or \code{"all"}.
#' @param available Available score column names.
#'
#' @return Filtered \code{edges}.
#' @noRd
.filter_edges_by_scores <- function(
  edges,
  score_min,
  score_max,
  score_combine,
  available
) {
  score_min <- .validate_named_score_bounds(score_min, "score_min", available)
  score_max <- .validate_named_score_bounds(score_max, "score_max", available)
  if (is.null(score_min) && is.null(score_max)) {
    return(edges)
  }

  cols <- union(names(score_min), names(score_max))
  shared <- intersect(names(score_min), names(score_max))
  for (col in shared) {
    if (score_min[[col]] >= score_max[[col]]) {
      cli::cli_abort(c(
        "x" = "Invalid score range for {.val {col}}.",
        "i" = "{.arg score_min} ({score_min[[col]]}) must be < {.arg score_max} ({score_max[[col]]})."
      ))
    }
  }

  preds <- lapply(cols, function(col) {
    values <- edges[[col]]
    keep <- !is.na(values)
    if (!is.null(score_min) && col %in% names(score_min)) {
      keep <- keep & (values >= score_min[[col]])
    }
    if (!is.null(score_max) && col %in% names(score_max)) {
      keep <- keep & (values < score_max[[col]])
    }
    keep
  })
  keep_rows <- if (identical(score_combine, "any")) {
    Reduce(`|`, preds)
  } else {
    Reduce(`&`, preds)
  }
  return(edges[keep_rows, , drop = FALSE])
}

#' Default cache directory for interaction databases
#'
#' @return Character path under \code{tools::R_user_dir("pixelatorR", "cache")}.
#' @export
interaction_database_cache_dir <- function() {
  return(file.path(
    tools::R_user_dir("pixelatorR", which = "cache"),
    "interaction_databases"
  ))
}

#' Normalise an edge table to the canonical interaction-database schema
#'
#' @param edges Data frame with at least two UniProt columns.
#' @param a_col,b_col Column names for the two ends.
#' @param score_cols Optional character vector of score columns to retain.
#' @param evidence_col Optional evidence column.
#' @param resource Resource name.
#' @param resource_version Version string.
#' @return Tibble with canonical columns.
#' @export
normalise_interaction_edges <- function(
  edges,
  a_col = "uniprot_a",
  b_col = "uniprot_b",
  score_cols = NULL,
  evidence_col = NULL,
  resource,
  resource_version
) {
  assert_class(edges, c("data.frame", "tbl_df"))
  assert_single_value(a_col, type = "string")
  assert_single_value(b_col, type = "string")
  assert_vector(
    score_cols,
    type = "character",
    n = 1,
    allow_null = TRUE
  )
  assert_single_value(evidence_col, type = "string", allow_null = TRUE)
  assert_single_value(resource, type = "string")
  assert_single_value(resource_version, type = "string")
  assert_col_in_data(a_col, edges)
  assert_col_in_data(b_col, edges)
  if (!is.null(score_cols)) {
    missing_score_cols <- setdiff(score_cols, names(edges))
    if (length(missing_score_cols) > 0) {
      cli::cli_abort(
        "Missing score columns: {.val {missing_score_cols}}."
      )
    }
  }

  a <- trimws(as.character(edges[[a_col]]))
  b <- trimws(as.character(edges[[b_col]]))
  keep <- !is.na(a) & !is.na(b) & a != "" & b != ""
  a <- a[keep]
  b <- b[keep]
  output <- tibble(
    uniprot_a = pmin(a, b),
    uniprot_b = pmax(a, b)
  )
  if (!is.null(score_cols)) {
    scores <- lapply(edges[keep, score_cols, drop = FALSE], function(score) {
      # as.character first so factor scores are not silently converted to codes
      as.numeric(as.character(score))
    })
    output <- bind_cols(output, as_tibble(scores))
  }
  output <- bind_cols(output, tibble(
    evidence = if (!is.null(evidence_col) && evidence_col %in% names(edges)) {
      as.character(edges[[evidence_col]][keep])
    } else {
      NA_character_
    },
    resource = resource,
    resource_version = resource_version
  ))
  return(unique(output))
}

#' Write a slim interaction database RDS
#'
#' @param edges Canonical edge tibble.
#' @param database Database key.
#' @param version Version label used in the filename.
#' @param cache_dir Cache directory.
#' @param score_columns Character vector of numeric score column names present
#'   in \code{edges}. Use \code{NULL} or \code{character(0)} when the database
#'   has no scores.
#' @param source_url Optional provenance URL.
#' @param license Optional license string.
#' @param citation Optional citation string.
#' @return Path to the written RDS.
#' @seealso \code{\link{load_interaction_database}},
#'   \code{\link{build_all_interaction_databases}}
#' @export
save_interaction_database <- function(
  edges,
  database = c("string", "biogrid", "corum", "omnipath", "alphafold"),
  version,
  cache_dir = interaction_database_cache_dir(),
  score_columns = NULL,
  source_url = NULL,
  license = NULL,
  citation = NULL
) {
  database <- match.arg(database)
  assert_class(edges, c("data.frame", "tbl_df"))
  assert_single_value(version, type = "string")
  assert_single_value(cache_dir, type = "string")
  assert_single_value(source_url, type = "string", allow_null = TRUE)
  assert_single_value(license, type = "string", allow_null = TRUE)
  assert_single_value(citation, type = "string", allow_null = TRUE)

  required_cols <- c(
    "uniprot_a", "uniprot_b", "evidence", "resource", "resource_version"
  )
  missing_cols <- setdiff(required_cols, names(edges))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database edges.",
      "i" = "Missing required columns: {.val {missing_cols}}"
    ))
  }
  if (is.null(score_columns)) {
    score_columns <- character()
  } else if (length(score_columns) > 0) {
    assert_vector(score_columns, type = "character", n = 1)
  } else if (!is.character(score_columns)) {
    cli::cli_abort("{.arg score_columns} must be a character vector or NULL.")
  }
  if (anyDuplicated(score_columns)) {
    cli::cli_abort("{.arg score_columns} must not contain duplicates.")
  }
  reserved_scores <- intersect(score_columns, required_cols)
  if (length(reserved_scores) > 0) {
    cli::cli_abort(c(
      "x" = "{.arg score_columns} collides with required columns.",
      "i" = "Invalid names: {.val {reserved_scores}}"
    ))
  }
  missing_score_cols <- setdiff(score_columns, names(edges))
  if (length(missing_score_cols) > 0) {
    cli::cli_abort(
      "Missing score columns: {.val {missing_score_cols}}."
    )
  }
  if (length(score_columns) > 0) {
    non_numeric_scores <- score_columns[
      !vapply(edges[score_columns], is.numeric, logical(1))
    ]
    if (length(non_numeric_scores) > 0) {
      cli::cli_abort(
        "Score columns must be numeric: {.val {non_numeric_scores}}."
      )
    }
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- .interaction_database_rds_path(database, version, cache_dir)
  meta <- list(
    database = database,
    version = version,
    n_edges = nrow(edges),
    score_columns = score_columns,
    source_url = source_url,
    license = license,
    citation = citation,
    built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    file = basename(path)
  )
  saveRDS(list(edges = edges, meta = meta), path)
  latest <- .interaction_database_rds_path(database, "latest", cache_dir)
  if (!identical(path, latest)) {
    file.copy(path, latest, overwrite = TRUE)
  }
  return(path)
}

#' Load a slim interaction database
#'
#' @param database One of \code{"string"}, \code{"biogrid"}, \code{"corum"},
#'   \code{"omnipath"}, or \code{"alphafold"}.
#' @param version Version label or \code{"latest"}.
#' @param cache_dir Cache directory.
#' @return List with \code{edges} (tibble) and \code{meta}.
#' @seealso \code{\link{save_interaction_database}},
#'   \code{\link{extract_panel_interactions}}
#' @export
load_interaction_database <- function(
  database = c("string", "biogrid", "corum", "omnipath", "alphafold"),
  version = "latest",
  cache_dir = interaction_database_cache_dir()
) {
  database <- match.arg(database)
  assert_single_value(version, type = "string")
  assert_single_value(cache_dir, type = "string")

  path <- .interaction_database_rds_path(database, version, cache_dir)
  assert_file_exists(path)
  obj <- readRDS(path)
  assert_class(obj, "list")
  if (is.null(obj$edges)) {
    cli::cli_abort("Invalid interaction database object at {.path {path}}.")
  }
  assert_class(obj$edges, c("data.frame", "tbl_df"))
  required_cols <- c(
    "uniprot_a", "uniprot_b", "evidence", "resource", "resource_version"
  )
  missing_cols <- setdiff(required_cols, names(obj$edges))
  if (length(missing_cols) > 0) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database edges at {.path {path}}.",
      "i" = "Missing required columns: {.val {missing_cols}}"
    ))
  }
  score_columns <- obj$meta$score_columns %||% character()
  if (!is.character(score_columns)) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database meta at {.path {path}}.",
      "i" = "{.field score_columns} must be a character vector."
    ))
  }
  missing_score_cols <- setdiff(score_columns, names(obj$edges))
  if (length(missing_score_cols) > 0) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database edges at {.path {path}}.",
      "i" = "Missing declared score columns: {.val {missing_score_cols}}"
    ))
  }
  if (length(score_columns) > 0) {
    non_numeric_scores <- score_columns[
      !vapply(obj$edges[score_columns], is.numeric, logical(1))
    ]
    if (length(non_numeric_scores) > 0) {
      cli::cli_abort(c(
        "x" = "Invalid interaction database edges at {.path {path}}.",
        "i" = "Score columns must be numeric: {.val {non_numeric_scores}}"
      ))
    }
  }
  obj$meta$score_columns <- score_columns
  return(obj)
}

#' Create a marker-UniProt map from a Seurat / PNA object
#'
#' Reads feature metadata from a PNA (or other) assay and expands
#' semicolon-separated UniProt accessions to long form for use with
#' \code{\link{extract_panel_interactions}}.
#'
#' Equivalent to:
#' \preformatted{
#' object[[assay]][[]] %>%
#'   as_tibble(rownames = "marker") %>%
#'   select(marker, uniprot_id) %>%
#'   separate_rows(uniprot_id, sep = ";")
#' }
#'
#' @param object A \code{Seurat} object with feature metadata on the chosen assay.
#' @param assay Assay name. Defaults to \code{DefaultAssay(object)}.
#' @param uniprot_col Feature-metadata column with UniProt accessions
#'   (semicolon-separated values are expanded).
#' @return A tibble with columns \code{marker} and \code{uniprot_id}.
#'
#' @seealso \code{\link{extract_panel_interactions}}
#'
#' @examples
#' \dontrun{
#' marker_map <- create_marker_uniprot_map(pg_data, assay = "PNA")
#' }
#' @export
create_marker_uniprot_map <- function(
  object,
  assay = NULL,
  uniprot_col = "uniprot_id"
) {
  assert_class(object, "Seurat")
  assert_single_value(assay, type = "string", allow_null = TRUE)
  assert_single_value(uniprot_col, type = "string")

  assay <- assay %||% DefaultAssay(object)
  if (!assay %in% Assays(object)) {
    cli::cli_abort("Assay {.val {assay}} not found in {.cls Seurat} object.")
  }

  feat_meta <- object[[assay]][[]]
  if (!uniprot_col %in% colnames(feat_meta)) {
    cli::cli_abort(c(
      "x" = "Feature metadata for assay {.val {assay}} lacks column {.val {uniprot_col}}."
    ))
  }

  return(
    as_tibble(feat_meta, rownames = "marker") %>%
      select(marker, uniprot_id = all_of(uniprot_col)) %>%
      mutate(
        marker = as.character(marker),
        uniprot_id = as.character(uniprot_id)
      ) %>%
      separate_rows(uniprot_id, sep = ";") %>%
      mutate(uniprot_id = trimws(uniprot_id)) %>%
      filter(!is.na(uniprot_id), uniprot_id != "") %>%
      distinct(marker, uniprot_id)
  )
}

#' Extract known database interactions for a marker panel
#'
#' Maps panel markers to UniProt accessions, loads a slim interaction database,
#' and returns undirected marker pairs with a registered database entry after
#' optional score / network filters.
#' UniProt self-interactions are returned only as marker self-pairs. When
#' multiple marker names map to the same UniProt accession, a self-interaction
#' does not create interactions between those different marker names.
#' Conversely, interactions between different UniProt accessions are not
#' returned as marker self-pairs. Marker names are ordered alphabetically, and
#' the corresponding UniProt accessions are reordered with them.
#'
#' @param markers Character vector of panel marker names.
#' @param database Interaction database key.
#' @param marker_uniprot_map Data frame/tibble with columns \code{marker} and
#'   \code{uniprot_id} (long form; multiple rows per marker allowed). Build with
#'   \code{\link{create_marker_uniprot_map}}.
#' @param score_min Named numeric vector of inclusive lower bounds
#'   (\code{column >= value}), or \code{NULL} to skip. Names must match score
#'   columns in the loaded database (see \code{meta$score_columns}).
#' @param score_max Named numeric vector of exclusive upper bounds
#'   (\code{column < value}), or \code{NULL} to skip.
#' @param score_combine How to combine predicates across score columns:
#'   \code{"any"} (OR, default) or \code{"all"} (AND). Bounds on the same
#'   column are always ANDed into one per-column predicate.
#' @param string_network For STRING: \code{"physical"} or \code{"full"}.
#' @param cache_dir Interaction database cache directory.
#' @param version Database version label (\code{"latest"} or a built version).
#' @return A tibble of panel edges with \code{marker_1}, \code{marker_2},
#'   \code{uniprot_a}, \code{uniprot_b}, \code{in_db = TRUE}, native score
#'   columns from the database, plus \code{evidence}, \code{resource}, and
#'   \code{resource_version}. Pass \code{marker_1}/\code{marker_2} columns to
#'   \code{ColocalizationHeatmap(highlight_pairs = ...)}.
#'
#' @section Score filtering:
#' Slim caches keep native score columns rather than a single synthetic
#' \code{score}. Discover them with
#' \code{load_interaction_database(...)$meta$score_columns}:
#'
#' | Database | Score columns |
#' | --- | --- |
#' | STRING | \code{combined_score} (classic 0-1000 scale) |
#' | AlphaFold | \code{ipSAE}, \code{pDockQ2} (typically 0-1) |
#' | BioGRID, CORUM, OmniPath | none |
#'
#' - \code{score_min}: keep rows where \code{column >= value}
#' - \code{score_max}: keep rows where \code{column < value}
#' - Same column in both: \code{(x >= min) & (x < max)}
#' - \code{score_combine}: how predicates **across columns** combine
#'   (\code{"any"} = OR, \code{"all"} = AND)
#' - Both args must be fully named numeric vectors when non-\code{NULL}
#' - Non-\code{NULL} thresholds on a database with no score columns error
#'
#' Examples:
#'
#' ```r
#' # STRING medium confidence
#' extract_panel_interactions(
#'   markers, "string", marker_map,
#'   score_min = c(combined_score = 400)
#' )
#'
#' # AlphaFold: keep if ipSAE >= 0.6 OR pDockQ2 >= 0.23
#' extract_panel_interactions(
#'   markers, "alphafold", marker_map,
#'   score_min = c(ipSAE = 0.6, pDockQ2 = 0.23),
#'   score_combine = "any"
#' )
#'
#' # AlphaFold: both minima must pass
#' extract_panel_interactions(
#'   markers, "alphafold", marker_map,
#'   score_min = c(ipSAE = 0.6, pDockQ2 = 0.23),
#'   score_combine = "all"
#' )
#'
#' # Range on ipSAE AND min on pDockQ2:
#' # (ipSAE >= 0.6 & ipSAE < 0.99) & (pDockQ2 >= 0.23)
#' extract_panel_interactions(
#'   markers, "alphafold", marker_map,
#'   score_min = c(ipSAE = 0.6, pDockQ2 = 0.23),
#'   score_max = c(ipSAE = 0.99),
#'   score_combine = "all"
#' )
#' ```
#'
#' @seealso \code{\link{ColocalizationHeatmap}},
#'   \code{\link{create_marker_uniprot_map}},
#'   \code{\link{load_interaction_database}},
#'   \code{\link{build_all_interaction_databases}}
#'
#' @examples
#' \dontrun{
#' marker_map <- create_marker_uniprot_map(pg_data, assay = "PNA")
#' highlight_pairs <- extract_panel_interactions(
#'   markers = c("CD3e", "CD4", "CD8"),
#'   database = "string",
#'   marker_uniprot_map = marker_map,
#'   score_min = c(combined_score = 400)
#' ) %>%
#'   select(marker_1, marker_2)
#'
#' ColocalizationHeatmap(
#'   coloc_summary,
#'   type = "dots",
#'   highlight_pairs = highlight_pairs
#' )
#' }
#' @export
extract_panel_interactions <- function(
  markers,
  database = c("string", "biogrid", "corum", "omnipath", "alphafold"),
  marker_uniprot_map,
  score_min = NULL,
  score_max = NULL,
  score_combine = c("any", "all"),
  string_network = c("physical", "full"),
  cache_dir = interaction_database_cache_dir(),
  version = "latest"
) {
  database <- match.arg(database)
  score_combine <- match.arg(score_combine)
  string_network <- match.arg(string_network)
  assert_vector(markers, type = "character", n = 1)
  assert_class(marker_uniprot_map, c("data.frame", "tbl_df"))
  assert_col_in_data("marker", marker_uniprot_map)
  assert_col_in_data("uniprot_id", marker_uniprot_map)
  assert_col_class("marker", marker_uniprot_map, classes = "character")
  assert_col_class("uniprot_id", marker_uniprot_map, classes = "character")
  assert_single_value(cache_dir, type = "string")
  assert_single_value(version, type = "string")

  db <- load_interaction_database(
    database = database,
    version = version,
    cache_dir = cache_dir
  )
  edges <- db$edges
  score_cols <- db$meta$score_columns %||% character()

  # Zero-row result matching the loaded database score schema.
  empty_panel_interactions <- function() {
    score_tbl <- as_tibble(
      setNames(rep(list(numeric()), length(score_cols)), score_cols)
    )
    return(bind_cols(
      tibble(
        marker_1 = character(),
        marker_2 = character(),
        uniprot_a = character(),
        uniprot_b = character(),
        in_db = logical()
      ),
      score_tbl,
      tibble(
        evidence = character(),
        resource = character(),
        resource_version = character()
      )
    ))
  }

  markers <- unique(markers)
  map <- marker_uniprot_map %>%
    filter(
      marker %in% markers,
      !is.na(uniprot_id),
      uniprot_id != ""
    ) %>%
    select(marker, uniprot_id) %>%
    mutate(uniprot_id = trimws(uniprot_id)) %>%
    distinct(marker, uniprot_id)
  if (nrow(map) == 0) {
    return(empty_panel_interactions())
  }

  if (database == "string") {
    edges <- edges %>%
      filter(evidence == string_network)
  }
  edges <- .filter_edges_by_scores(
    edges = edges,
    score_min = score_min,
    score_max = score_max,
    score_combine = score_combine,
    available = score_cols
  )

  panel_uniprot <- unique(map$uniprot_id)
  edges <- edges %>%
    filter(
      uniprot_a %in% panel_uniprot,
      uniprot_b %in% panel_uniprot
    )

  if (nrow(edges) == 0) {
    return(empty_panel_interactions())
  }

  map_a <- map %>% rename(marker_1 = marker)
  map_b <- map %>% rename(marker_2 = marker)
  out <- edges %>%
    inner_join(map_a, by = c("uniprot_a" = "uniprot_id")) %>%
    inner_join(map_b, by = c("uniprot_b" = "uniprot_id")) %>%
    # Homodimers -> marker self-pairs only; hetero edges -> distinct markers only
    filter(
      (uniprot_a == uniprot_b & marker_1 == marker_2) |
        (uniprot_a != uniprot_b & marker_1 != marker_2)
    ) %>%
    mutate(
      swap = marker_1 > marker_2,
      marker_1_ordered = if_else(swap, marker_2, marker_1),
      marker_2_ordered = if_else(swap, marker_1, marker_2),
      uniprot_a_ordered = if_else(swap, uniprot_b, uniprot_a),
      uniprot_b_ordered = if_else(swap, uniprot_a, uniprot_b)
    ) %>%
    select(
      marker_1 = marker_1_ordered,
      marker_2 = marker_2_ordered,
      uniprot_a = uniprot_a_ordered,
      uniprot_b = uniprot_b_ordered,
      all_of(score_cols),
      evidence,
      resource,
      resource_version
    ) %>%
    group_by(marker_1, marker_2) %>%
    summarise(
      uniprot_a = uniprot_a[[1]],
      uniprot_b = uniprot_b[[1]],
      in_db = TRUE,
      across(
        all_of(score_cols),
        ~ suppressWarnings(max(.x, na.rm = TRUE))
      ),
      evidence = evidence[[1]],
      resource = resource[[1]],
      resource_version = resource_version[[1]],
      .groups = "drop"
    )
  if (length(score_cols) > 0) {
    out <- out %>%
      mutate(across(
        all_of(score_cols),
        ~ if_else(is.infinite(.x), NA_real_, .x)
      ))
  }
  return(out)
}

#' Build STRING physical / full link database (human)
#'
#' Maintainer helper: reads STRING release files (downloading them into
#' \code{raw_dir} when missing) and writes a slim UniProt edge cache via
#' \code{\link{save_interaction_database}}.
#'
#' @param raw_dir Directory for STRING download files (links + aliases).
#' @param version STRING version label (default \code{"12.0"}).
#' @param cache_dir Output interaction database cache.
#' @param species NCBI taxon (default human 9606).
#' @param include_full If \code{TRUE}, also include the full STRING network.
#' @return Path to saved RDS.
#' @export
build_string_database <- function(
  raw_dir = file.path(.default_raw_cache(), "string"),
  version = "12.0",
  cache_dir = interaction_database_cache_dir(),
  species = 9606,
  include_full = FALSE
) {
  rlang::check_installed("data.table")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  string_download_url <- function(fname, network = c("aliases", "physical", "full")) {
    network <- match.arg(network)
    prefix <- switch(network,
      aliases = paste0("protein.aliases.v", version),
      physical = paste0("protein.physical.links.v", version),
      full = paste0("protein.links.v", version)
    )
    paste0("https://stringdb-downloads.org/download/", prefix, "/", fname)
  }

  # Download a STRING release file into raw_dir when missing.
  ensure_string_file <- function(fname, network) {
    path <- file.path(raw_dir, fname)
    .download_if_missing(string_download_url(fname, network), path)
  }

  aliases_fname <- sprintf("%s.protein.aliases.v%s.txt.gz", species, version)
  physical_fname <- sprintf(
    "%s.protein.physical.links.v%s.txt.gz", species, version
  )
  full_fname <- sprintf("%s.protein.links.v%s.txt.gz", species, version)

  aliases_gz <- ensure_string_file(aliases_fname, "aliases")
  physical_gz <- ensure_string_file(physical_fname, "physical")
  full_gz <- if (isTRUE(include_full)) {
    ensure_string_file(full_fname, "full")
  } else {
    NULL
  }

  aliases <- data.table::fread(
    aliases_gz,
    header = TRUE,
    sep = "\t",
    quote = "",
    showProgress = FALSE
  )
  data.table::setnames(aliases, c("string_protein_id", "alias", "source"))
  # One UniProt per STRING protein: prefer the primary Swiss-Prot accession.
  # UniProt_AC alone includes secondaries (often listed first, e.g. O43746 before
  # P20701 for ITGAL). Expanding all ACs is correct but ~20x slower; prefer
  # Ensembl_HGNC_uniprot_ids / neXtProt primary, then classic Swiss-Prot AC.
  ac <- aliases[
    source == "UniProt_AC" & !is.na(alias) & alias != "",
    .(string_protein_id, uniprot = alias)
  ]
  preferred <- unique(rbind(
    aliases[
      source == "Ensembl_HGNC_uniprot_ids" & !is.na(alias) & alias != "",
      .(string_protein_id, uniprot = alias)
    ],
    aliases[
      source == "UniProt_DR_neXtProt" & grepl("^NX_", alias),
      .(string_protein_id, uniprot = sub("^NX_", "", alias))
    ]
  ))
  primary <- unique(merge(
    ac,
    preferred,
    by = c("string_protein_id", "uniprot")
  ))
  primary <- primary[, .(uniprot = uniprot[[1]]), by = string_protein_id]
  fallback_ids <- setdiff(unique(ac$string_protein_id), primary$string_protein_id)
  if (length(fallback_ids) > 0) {
    fallback <- ac[string_protein_id %in% fallback_ids][
      ,
      {
        swiss <- uniprot[grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", uniprot)]
        .(uniprot = if (length(swiss)) swiss[[1]] else uniprot[[1]])
      },
      by = string_protein_id
    ]
    up <- rbind(primary, fallback)
  } else {
    up <- primary
  }

  read_links <- function(path, network_label) {
    if (length(path) != 1 || is.na(path) || !file.exists(path)) {
      return(NULL)
    }
    links <- data.table::fread(path, header = TRUE, sep = " ", showProgress = FALSE)
    data.table::setnames(links, tolower(gsub("#", "", names(links))))
    if (!all(c("protein1", "protein2", "combined_score") %in% names(links))) {
      cli::cli_abort("Unexpected STRING links columns in {.path {path}}.")
    }
    m2 <- as_tibble(links) %>%
      inner_join(as_tibble(up), by = c("protein1" = "string_protein_id")) %>%
      rename(uniprot_a = uniprot) %>%
      inner_join(as_tibble(up), by = c("protein2" = "string_protein_id")) %>%
      rename(uniprot_b = uniprot)
    edges <- normalise_interaction_edges(
      m2,
      a_col = "uniprot_a",
      b_col = "uniprot_b",
      score_cols = "combined_score",
      evidence_col = NULL,
      resource = paste0("string_", network_label),
      resource_version = version
    )
    edges$evidence <- network_label
    edges
  }

  phys <- read_links(physical_gz, "physical")
  full <- if (isTRUE(include_full)) read_links(full_gz, "full") else NULL
  if (is.null(phys) && is.null(full)) {
    cli::cli_abort(
      "No STRING link files found (expected physical and/or full links .txt.gz)."
    )
  }
  edges <- bind_rows(phys, full)

  return(save_interaction_database(
    edges = edges,
    database = "string",
    version = version,
    cache_dir = cache_dir,
    score_columns = "combined_score",
    source_url = "https://stringdb-downloads.org/download/",
    license = paste(
      "Creative Commons Attribution 4.0 International (CC BY 4.0);",
      "attribution required"
    ),
    citation = paste(
      "Szklarczyk D, Kirsch R, Koutrouli M, et al. (2023).",
      "The STRING database in 2023: protein-protein association networks",
      "and functional enrichment analyses for any sequenced genome of",
      "interest. Nucleic Acids Research 51(D1):D638-D646.",
      "https://doi.org/10.1093/nar/gkac1000"
    )
  ))
}

#' Build BioGRID physical interaction database (human)
#'
#' Downloads \code{BIOGRID-MV-Physical-LATEST.tab3.zip} into the BioGRID raw
#' cache when \code{raw_file} is missing.
#'
#' @param raw_file Path to BioGRID MV-Physical or organism tab3 file.
#'   If \code{NULL} or missing on disk, the latest release is downloaded.
#' @param version BioGRID release label (default \code{"latest"} when downloading).
#' @param cache_dir Output interaction database cache.
#' @return Path to saved RDS.
#' @export
build_biogrid_database <- function(
  raw_file = NULL,
  version = "latest",
  cache_dir = interaction_database_cache_dir()
) {
  rlang::check_installed("data.table")
  raw_dir <- file.path(.default_raw_cache(), "biogrid")

  # Download a zip and return the first extracted file matching pattern.
  download_zip_extract <- function(url, dest_dir, pattern) {
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    existing <- list.files(dest_dir, pattern = pattern, full.names = TRUE)
    if (length(existing) > 0) {
      return(existing[[1]])
    }
    zip_path <- file.path(dest_dir, basename(sub("[?].*$", "", url)))
    if (!grepl("[.]zip$", zip_path, ignore.case = TRUE)) {
      zip_path <- paste0(zip_path, ".zip")
    }
    .download_if_missing(url, zip_path)
    utils::unzip(zip_path, exdir = dest_dir)
    extracted <- list.files(dest_dir, pattern = pattern, full.names = TRUE)
    if (length(extracted) < 1) {
      cli::cli_abort(
        "No file matching {.val {pattern}} found after unzipping {.path {zip_path}}."
      )
    }
    return(extracted[[1]])
  }

  if (is.null(raw_file)) {
    existing <- list.files(raw_dir, pattern = "[.]tab3[.]txt$", full.names = TRUE)
    raw_file <- if (length(existing) > 0) existing[[1]] else NA_character_
  }
  if (length(raw_file) != 1 || is.na(raw_file) || !file.exists(raw_file)) {
    raw_file <- download_zip_extract(
      url = paste0(
        "https://downloads.thebiogrid.org/Download/BioGRID/",
        "Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
      ),
      dest_dir = raw_dir,
      pattern = "[.]tab3[.]txt$"
    )
  }
  dt <- data.table::fread(raw_file, sep = "\t", quote = "", showProgress = FALSE)
  a_col <- .first_present_col(
    dt,
    c(
      "SWISS-PROT.*Interactor A",
      "SWISS-PROT Accessions Interactor A"
    ),
    regex = TRUE
  )
  b_col <- .first_present_col(
    dt,
    c(
      "SWISS-PROT.*Interactor B",
      "SWISS-PROT Accessions Interactor B"
    ),
    regex = TRUE
  )
  org_a <- .first_present_col(dt, "Organism ID Interactor A", regex = TRUE)
  org_b <- .first_present_col(dt, "Organism ID Interactor B", regex = TRUE)
  exp_col <- .first_present_col(dt, "^Experimental System$", regex = TRUE)
  if (is.na(a_col) || is.na(b_col)) {
    cli::cli_abort("Could not find SWISS-PROT columns in {.path {raw_file}}.")
  }

  if (!is.na(org_a) && !is.na(org_b)) {
    dt <- dt[as.character(dt[[org_a]]) == "9606" & as.character(dt[[org_b]]) == "9606"]
  }

  first_ac <- function(x) {
    v <- trimws(sub("\\|.*$", "", as.character(x)))
    v[v == "-" | v == ""] <- NA_character_
    v
  }
  evidence <- if (!is.na(exp_col)) as.character(dt[[exp_col]]) else NA_character_
  pairs <- data.frame(
    uniprot_a = first_ac(dt[[a_col]]),
    uniprot_b = first_ac(dt[[b_col]]),
    evidence = evidence,
    stringsAsFactors = FALSE
  )
  pairs <- pairs[!is.na(pairs$uniprot_a) & !is.na(pairs$uniprot_b), , drop = FALSE]
  edges <- normalise_interaction_edges(
    pairs,
    a_col = "uniprot_a",
    b_col = "uniprot_b",
    evidence_col = "evidence",
    resource = "biogrid",
    resource_version = version
  )

  return(save_interaction_database(
    edges = edges,
    database = "biogrid",
    version = version,
    cache_dir = cache_dir,
    score_columns = NULL,
    source_url = "https://downloads.thebiogrid.org/BioGRID/Latest-Release/",
    license = paste(
      "MIT License; retain the copyright and permission notices from",
      "BioGRID's LICENSE.txt"
    ),
    citation = paste(
      "Stark C, Breitkreutz BJ, Reguly T, Boucher L, Breitkreutz A,",
      "Tyers M. (2006). BioGRID: a general repository for interaction",
      "datasets. Nucleic Acids Research 34(Database issue):D535-D539.",
      "https://doi.org/10.1093/nar/gkj109; also cite original",
      "contributing publications where applicable"
    )
  ))
}

#' Build CORUM co-membership database (human)
#'
#' Downloads OmniPath-served CORUM complexes when \code{corum_file} is missing.
#' Prefer an official CORUM file when present. Complexes with a single UniProt
#' accession become homomer edges (\code{uniprot_a == uniprot_b}).
#'
#' @param corum_file Path to OmniPath CORUM complexes TSV or official coreComplexes.
#' @param version Version label.
#' @param cache_dir Output interaction database cache.
#' @return Path to saved RDS.
#' @export
build_corum_database <- function(
  corum_file = file.path(
    .default_raw_cache(), "corum", "omnipath_complexes_corum.tsv"
  ),
  version = "omnipath_corum",
  cache_dir = interaction_database_cache_dir()
) {
  rlang::check_installed("data.table")
  if (!file.exists(corum_file)) {
    .download_if_missing(
      "https://omnipathdb.org/complexes?databases=CORUM",
      corum_file
    )
  }
  first <- readLines(corum_file, n = 5, warn = FALSE)
  if (any(grepl("<!DOCTYPE html|<html", first, ignore.case = TRUE))) {
    cli::cli_abort("CORUM file looks like HTML, not data: {.path {corum_file}}.")
  }

  dt <- data.table::fread(corum_file, sep = "\t", quote = "", showProgress = FALSE)
  comp_col <- .first_present_col(
    dt,
    c(
      "components",
      "subunits_uniprot_id",
      "subunits",
      "subunits(UniProt IDs)"
    )
  )
  if (is.na(comp_col)) {
    cli::cli_abort("No components column in {.path {corum_file}}.")
  }
  name_col <- .first_present_col(
    dt,
    c("name", "complex_name", "ComplexName")
  )

  pairs <- bind_rows(
    tibble(
      uniprot_a = character(),
      uniprot_b = character(),
      evidence = character()
    ),
    lapply(seq_len(nrow(dt)), function(i) {
      ids <- unique(trimws(unlist(strsplit(
        as.character(dt[[comp_col]][[i]]),
        "[_;,]"
      ))))
      ids <- ids[ids != "" & !is.na(ids)]
      evid <- if (!is.na(name_col)) {
        as.character(dt[[name_col]][[i]])
      } else {
        NA_character_
      }
      if (length(ids) < 1) {
        return(NULL)
      }
      # Single-accession complexes (homomers) become U-U edges
      if (length(ids) == 1) {
        return(tibble(
          uniprot_a = ids,
          uniprot_b = ids,
          evidence = evid
        ))
      }
      grid <- utils::combn(ids, 2)
      tibble(
        uniprot_a = grid[1, ],
        uniprot_b = grid[2, ],
        evidence = evid
      )
    })
  )
  edges <- normalise_interaction_edges(
    pairs,
    evidence_col = "evidence",
    resource = "corum",
    resource_version = version
  )

  return(save_interaction_database(
    edges = edges,
    database = "corum",
    version = version,
    cache_dir = cache_dir,
    score_columns = NULL,
    source_url = paste(
      "https://omnipathdb.org/complexes?databases=CORUM;",
      "https://mips.helmholtz-muenchen.de/corum/"
    ),
    license = paste(
      "Creative Commons Attribution 4.0 International (CC BY 4.0);",
      "attribution required"
    ),
    citation = paste(
      "Steinkamp R, Tsitsiridis G, Brauner B, Montrone C, Fobo G,",
      "Frishman G, Avram S, Oprea TI, Ruepp A. (2025). CORUM in 2024:",
      "protein complexes as drug targets. Nucleic Acids Research",
      "53(D1):D651-D657. https://doi.org/10.1093/nar/gkae1033"
    )
  ))
}

#' Build OmniPath interaction database
#'
#' @param interactions_file Path to OmniPath interactions TSV.
#' @param version Version label.
#' @param cache_dir Output interaction database cache.
#' @param license OmniPath license filter used when the file was downloaded
#'   (\code{"commercial"} recommended for product builds).
#' @return Path to saved RDS.
#' @export
build_omnipath_database <- function(
  interactions_file = NULL,
  version = "webservice",
  cache_dir = interaction_database_cache_dir(),
  license = c("commercial", "academic", "unknown")
) {
  rlang::check_installed("data.table")
  license <- match.arg(license)
  if (is.null(interactions_file)) {
    interactions_file <- file.path(
      .default_raw_cache(), "omnipath",
      paste0("interactions_omnipath_", license, ".tsv")
    )
    fallback <- file.path(
      .default_raw_cache(), "omnipath", "interactions_omnipath.tsv"
    )
    if (
      license == "unknown" &&
        !file.exists(interactions_file) &&
        file.exists(fallback)
    ) {
      interactions_file <- fallback
    }
  }
  if (!file.exists(interactions_file)) {
    if (license == "unknown") {
      cli::cli_abort(
        "Missing OmniPath interactions file: {.path {interactions_file}}."
      )
    }
    url <- paste0(
      "https://omnipathdb.org/interactions?",
      "datasets=omnipath&fields=sources,references,curation_effort&genesymbols=1",
      "&license=", license
    )
    .download_if_missing(url, interactions_file)
  }
  dt <- data.table::fread(
    interactions_file,
    sep = "\t",
    quote = "",
    showProgress = FALSE
  )
  src <- .first_present_col(dt, c("source", "uniprot_a"))
  tgt <- .first_present_col(dt, c("target", "uniprot_b"))
  if (is.na(src) || is.na(tgt)) {
    cli::cli_abort("Could not find source/target UniProt columns.")
  }
  evid_col <- intersect(c("sources", "references"), colnames(dt))
  evid <- if (length(evid_col)) {
    apply(dt[, evid_col, with = FALSE], 1, function(x) paste(x, collapse = "|"))
  } else {
    NA_character_
  }
  pairs <- data.frame(
    uniprot_a = as.character(dt[[src]]),
    uniprot_b = as.character(dt[[tgt]]),
    evidence = evid,
    stringsAsFactors = FALSE
  )
  edges <- normalise_interaction_edges(
    pairs,
    evidence_col = "evidence",
    resource = "omnipath",
    resource_version = paste0(version, "_", license)
  )

  return(save_interaction_database(
    edges = edges,
    database = "omnipath",
    version = paste0(version, "_", license),
    cache_dir = cache_dir,
    score_columns = NULL,
    source_url = "https://omnipathdb.org/interactions",
    license = paste0(
      "No single OmniPath data license; contributing-resource licenses ",
      "apply; web-service license filter=", license,
      ". Cite and comply with each contributing resource"
    ),
    citation = paste(
      "T\u00fcrei D, Valdeolivas A, Gul L, et al. (2021). Integrated intra-",
      "and intercellular signaling knowledge for multicellular omics",
      "analysis. Molecular Systems Biology 17:e9923.",
      "https://doi.org/10.15252/msb.20209923; also cite the contributing",
      "resources identified in the evidence/sources fields"
    )
  ))
}

#' Build AlphaFold DB high-confidence complex database
#'
#' Combines heterodimer and homodimer predictions. When no local CSVs are
#' present, downloads both NVIDIA/AFDB metadata tables from EBI FTP into
#' \code{raw_dir} (heterodimer ~2 GB, homodimer ~6 GB).
#'
#' Homodimer metadata uses columns \code{uniprotAccession} and \code{taxId};
#' heterodimers use \code{uniprot_ac_1}/\code{uniprot_ac_2} and
#' \code{tax_id_1}/\code{tax_id_2}.
#'
#' @param heterodimer_file Optional FTP-derived heterodimer metadata CSV.
#'   If absent, uses panel API cache CSVs under \code{raw_dir}, then downloads
#'   the official heterodimer and homodimer metadata dumps if needed.
#' @param raw_dir Directory with AFDB cache files.
#' @param version Version label.
#' @param cache_dir Output interaction database cache.
#' @param ipsae_min Minimum ipSAE.
#' @param pdockq2_min Minimum pDockQ2.
#' @return Path to saved RDS.
#' @export
build_alphafold_database <- function(
  heterodimer_file = file.path(
    .default_raw_cache(), "alphafold_interactions",
    "afdb_heterodimers_human_panel_any.csv"
  ),
  raw_dir = file.path(.default_raw_cache(), "alphafold_interactions"),
  version = "afdb_nvda_highconf",
  cache_dir = interaction_database_cache_dir(),
  ipsae_min = 0.6,
  pdockq2_min = 0.23
) {
  # First matching column's values among exact name candidates, or NULL.
  coalesce_col <- function(df, candidates) {
    col <- .first_present_col(df, candidates, regex = FALSE)
    if (is.na(col)) {
      return(NULL)
    }
    return(df[[col]])
  }

  # Max of directional *_AB / *_BA scores, keeping a non-missing single side.
  max_directional_score <- function(x, y) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    ifelse(is.na(x), y, ifelse(is.na(y), x, pmax(x, y)))
  }

  parse_afdb_file <- function(f) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    # Homodimer NVIDIA schema: single UniProt + taxId
    if ("uniprotAccession" %in% names(df)) {
      if ("taxId" %in% names(df)) {
        df <- df[as.character(df$taxId) == "9606", , drop = FALSE]
      }
      a <- as.character(df$uniprotAccession)
      b <- a
    } else {
      if (all(c("tax_id_1", "tax_id_2") %in% names(df))) {
        df <- df[
          as.character(df$tax_id_1) == "9606" &
            as.character(df$tax_id_2) == "9606", ,
          drop = FALSE
        ]
      }
      a <- coalesce_col(df, c("uniprot_ac_1", "uniprot_a", "a"))
      b <- coalesce_col(df, c("uniprot_ac_2", "uniprot_b", "b"))
    }
    if (is.null(a) || is.null(b) || length(a) < 1) {
      return(NULL)
    }
    ipsae <- coalesce_col(df, c("ipSAE", "ipsae"))
    if (!is.null(ipsae)) {
      ipsae <- as.numeric(ipsae)
    } else if (all(c("ipSAE_AB", "ipSAE_BA") %in% names(df))) {
      ipsae <- max_directional_score(df$ipSAE_AB, df$ipSAE_BA)
    } else {
      ipsae <- rep(NA_real_, length(a))
    }
    pdq <- coalesce_col(df, c("pDockQ2", "pdockq2"))
    if (!is.null(pdq)) {
      pdq <- as.numeric(pdq)
    } else if (all(c("pDockQ2_AB", "pDockQ2_BA") %in% names(df))) {
      pdq <- max_directional_score(df$pDockQ2_AB, df$pDockQ2_BA)
    } else {
      pdq <- rep(NA_real_, length(a))
    }
    data.frame(
      uniprot_a = as.character(a),
      uniprot_b = as.character(b),
      ipSAE = ipsae,
      pDockQ2 = pdq,
      stringsAsFactors = FALSE
    )
  }

  hetero_dest <- file.path(raw_dir, "heterodimer_metadata.csv")
  homo_dest <- file.path(raw_dir, "homodimer_metadata.csv")
  nvda_base <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/alphafold/",
    "collaborations/nvda/"
  )

  files <- unique(c(
    heterodimer_file,
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    file.path(raw_dir, "afdb_heterodimers_panel_pairs.csv"),
    file.path(raw_dir, "afdb_homodimers_human_panel.csv"),
    hetero_dest,
    homo_dest
  ))
  files <- files[file.exists(files)]

  # Only download official dumps when nothing local is present. Do not fetch
  # a missing counterpart when a small panel CSV / fixture is already in use.
  if (length(files) < 1) {
    cli::cli_inform(c(
      "i" = paste(
        "No local AFDB CSVs found; downloading NVIDIA heterodimer (~2 GB)",
        "and homodimer (~6 GB) metadata."
      )
    ))
    .download_if_missing(
      paste0(nvda_base, "heterodimer_metadata.csv"),
      hetero_dest,
      min_timeout = 7200
    )
    .download_if_missing(
      paste0(nvda_base, "homodimer_metadata.csv"),
      homo_dest,
      min_timeout = 14400
    )
    files <- c(hetero_dest, homo_dest)
  }

  pieces <- lapply(files, parse_afdb_file)
  pairs <- bind_rows(pieces)
  if (nrow(pairs) < 1) {
    cli::cli_abort("No usable AFDB pairs found in {.path {raw_dir}}.")
  }
  keep <- (
    (!is.na(pairs$ipSAE) & pairs$ipSAE >= ipsae_min) |
      (!is.na(pairs$pDockQ2) & pairs$pDockQ2 >= pdockq2_min)
  )
  pairs <- pairs[keep, , drop = FALSE]
  edges <- normalise_interaction_edges(
    pairs,
    score_cols = c("ipSAE", "pDockQ2"),
    resource = "alphafold",
    resource_version = version
  )

  return(save_interaction_database(
    edges = edges,
    database = "alphafold",
    version = version,
    cache_dir = cache_dir,
    score_columns = c("ipSAE", "pDockQ2"),
    source_url = nvda_base,
    license = paste(
      "Creative Commons Attribution 4.0 International (CC BY 4.0);",
      "academic and commercial use permitted with attribution;",
      "EMBL-EBI Terms of Use also apply"
    ),
    citation = paste(
      "Bertoni D, Tsenkov MI, Maga\u00f1a P, et al. (2026). AlphaFold Protein",
      "Structure Database 2025: a redesigned interface and updated",
      "structural coverage. Nucleic Acids Research 54(D1):D358-D362.",
      "https://doi.org/10.1093/nar/gkaf1226; Han Y, Tsenkov MI,",
      "Venanzi NAE, et al. (2026). AlphaFold Database expands to",
      "proteome-scale quaternary structures. bioRxiv 2026.03.27.714458.",
      "https://doi.org/10.64898/2026.03.27.714458"
    )
  ))
}

#' Build all five interaction databases
#'
#' Maintainer helper that runs each \code{build_*_database()} writer. Missing
#' raw dumps are downloaded into the package raw cache. STRING is built with
#' both physical and full networks (\code{include_full = TRUE}). OmniPath
#' defaults to the commercial license filter.
#'
#' @param cache_dir Output interaction database cache.
#' @return Named list of RDS paths.
#' @seealso \code{\link{build_string_database}},
#'   \code{\link{build_biogrid_database}},
#'   \code{\link{build_corum_database}},
#'   \code{\link{build_omnipath_database}},
#'   \code{\link{build_alphafold_database}},
#'   \code{\link{extract_panel_interactions}}
#' @export
build_all_interaction_databases <- function(
  cache_dir = interaction_database_cache_dir()
) {
  return(list(
    string = build_string_database(
      cache_dir = cache_dir,
      include_full = TRUE
    ),
    biogrid = build_biogrid_database(cache_dir = cache_dir),
    corum = build_corum_database(cache_dir = cache_dir),
    omnipath = build_omnipath_database(
      cache_dir = cache_dir,
      license = "commercial"
    ),
    alphafold = build_alphafold_database(cache_dir = cache_dir)
  ))
}
