# Canonical edge schema: uniprot_a, uniprot_b (undirected, a <= b), optional
# native score columns and optional additional_columns. Provenance (resource,
# version, license, …) lives in meta only.
# Persistent cache is the edge RDS under interaction_database_cache_dir() only;
# raw vendor dumps are ephemeral staging during build_*_database().

#' Create an ephemeral staging directory for raw interaction-database dumps
#'
#' Used when builders download vendor files. Callers that pass their own
#' \code{raw_dir} / \code{raw_file} own those paths and they are not deleted.
#'
#' @param source Short source key (e.g. \code{"string"}).
#' @return Character path to a newly created empty directory under
#'   \code{tempdir()}.
#' @noRd
.interaction_database_staging_dir <- function(source) {
  assert_single_value(source, type = "string")
  staging <- tempfile(paste0("pixelatorR_", source, "_staging_"))
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  return(staging)
}

#' Finalize a raw staging directory after a build attempt
#'
#' On success with an owned staging dir, deletes it. On failure, leaves it and
#' informs so a retry can reuse the download. Never deletes caller-owned paths.
#'
#' @param staging_dir Directory path (may be \code{NULL}).
#' @param owned If \code{TRUE}, this package created \code{staging_dir}.
#' @param success If \code{TRUE}, the edge RDS was written successfully.
#' @return \code{NULL}, invisibly.
#' @noRd
.finalize_interaction_database_staging <- function(
  staging_dir,
  owned,
  success
) {
  if (!isTRUE(owned) || is.null(staging_dir) || !nzchar(staging_dir)) {
    return(invisible(NULL))
  }
  if (isTRUE(success)) {
    if (dir.exists(staging_dir)) {
      unlink(staging_dir, recursive = TRUE)
    }
    return(invisible(NULL))
  }
  if (dir.exists(staging_dir)) {
    cli::cli_inform(c(
      "i" = paste(
        "Leaving staging raw files at {.path {staging_dir}} for retry/debug.",
        "They are removed automatically after a successful rebuild."
      )
    ))
  }
  return(invisible(NULL))
}

#' Path to an interaction-database edge RDS file
#'
#' @param database Database key (e.g. \code{"string"}).
#' @param version Version label used in the filename (e.g. \code{"latest"}).
#' @param cache_dir Directory containing edge RDS caches.
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
#' Sole persistent store for built edge RDS files. Builders download raw vendor
#' dumps into ephemeral staging, write an RDS here, then delete the staging
#' files on success.
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
#' @param score_cols Optional character vector of numeric score columns to retain.
#' @param additional_cols Optional character vector of extra columns to copy
#'   through under their native names (non-score metadata such as
#'   \code{network} or \code{Experimental System}).
#' @return Tibble with canonical columns.
#' @export
normalise_interaction_edges <- function(
  edges,
  a_col = "uniprot_a",
  b_col = "uniprot_b",
  score_cols = NULL,
  additional_cols = NULL
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
  assert_vector(
    additional_cols,
    type = "character",
    n = 1,
    allow_null = TRUE
  )
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
  if (!is.null(additional_cols)) {
    if (anyDuplicated(additional_cols)) {
      cli::cli_abort("{.arg additional_cols} must not contain duplicates.")
    }
    missing_add_cols <- setdiff(additional_cols, names(edges))
    if (length(missing_add_cols) > 0) {
      cli::cli_abort(
        "Missing additional columns: {.val {missing_add_cols}}."
      )
    }
    overlap <- intersect(score_cols %||% character(), additional_cols)
    if (length(overlap) > 0) {
      cli::cli_abort(c(
        "x" = "{.arg additional_cols} overlaps {.arg score_cols}.",
        "i" = "Shared names: {.val {overlap}}"
      ))
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
  if (!is.null(score_cols) && length(score_cols) > 0) {
    scores <- lapply(edges[keep, score_cols, drop = FALSE], function(score) {
      # as.character first so factor scores are not silently converted to codes
      as.numeric(as.character(score))
    })
    output <- bind_cols(output, as_tibble(scores))
  }
  if (!is.null(additional_cols) && length(additional_cols) > 0) {
    extras <- edges[keep, additional_cols, drop = FALSE]
    # Character-coerce factors so unique() / joins stay stable
    extras[] <- lapply(extras, function(col) {
      if (is.factor(col) || is.character(col)) {
        as.character(col)
      } else {
        col
      }
    })
    output <- bind_cols(output, as_tibble(extras))
  }
  return(unique(output))
}

#' Write an interaction database edge RDS
#'
#' @param edges Canonical edge tibble.
#' @param database Database key.
#' @param version Version label used in the filename.
#' @param cache_dir Cache directory.
#' @param score_columns Character vector of numeric score column names present
#'   in \code{edges}. Use \code{NULL} or \code{character(0)} when the database
#'   has no scores.
#' @param additional_columns Character vector of extra non-score column names
#'   present in \code{edges}. Use \code{NULL} or \code{character(0)} when none.
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
  additional_columns = NULL,
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

  required_cols <- c("uniprot_a", "uniprot_b")
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
  if (is.null(additional_columns)) {
    additional_columns <- character()
  } else if (length(additional_columns) > 0) {
    assert_vector(additional_columns, type = "character", n = 1)
  } else if (!is.character(additional_columns)) {
    cli::cli_abort(
      "{.arg additional_columns} must be a character vector or NULL."
    )
  }
  if (anyDuplicated(additional_columns)) {
    cli::cli_abort("{.arg additional_columns} must not contain duplicates.")
  }
  reserved <- intersect(
    c(score_columns, additional_columns),
    required_cols
  )
  if (length(reserved) > 0) {
    cli::cli_abort(c(
      "x" = "Score/additional columns collide with required columns.",
      "i" = "Invalid names: {.val {reserved}}"
    ))
  }
  overlap <- intersect(score_columns, additional_columns)
  if (length(overlap) > 0) {
    cli::cli_abort(c(
      "x" = "{.arg additional_columns} overlaps {.arg score_columns}.",
      "i" = "Shared names: {.val {overlap}}"
    ))
  }
  missing_score_cols <- setdiff(score_columns, names(edges))
  if (length(missing_score_cols) > 0) {
    cli::cli_abort(
      "Missing score columns: {.val {missing_score_cols}}."
    )
  }
  missing_add_cols <- setdiff(additional_columns, names(edges))
  if (length(missing_add_cols) > 0) {
    cli::cli_abort(
      "Missing additional columns: {.val {missing_add_cols}}."
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
    resource = database,
    version = version,
    n_edges = nrow(edges),
    score_columns = score_columns,
    additional_columns = additional_columns,
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

#' Load an interaction database edge RDS
#'
#' @param database One of \code{"string"}, \code{"biogrid"}, \code{"corum"},
#'   \code{"omnipath"}, or \code{"alphafold"}.
#' @param version Version label or \code{"latest"}.
#' @param cache_dir Cache directory.
#' @param build_if_missing If \code{TRUE} and the edge RDS is missing, run the
#'   matching \code{build_*_database()} (download into ephemeral staging, write
#'   RDS, delete staging) then load. Defaults to \code{FALSE} so loads never
#'   surprise-download. Extra builder arguments can be passed via \code{...}.
#' @param ... Forwarded to \code{build_*_database()} when
#'   \code{build_if_missing = TRUE}.
#' @return List with \code{edges} (tibble) and \code{meta}.
#' @seealso \code{\link{save_interaction_database}},
#'   \code{\link{extract_panel_interactions}},
#'   \code{\link{build_all_interaction_databases}}
#' @export
load_interaction_database <- function(
  database = c("string", "biogrid", "corum", "omnipath", "alphafold"),
  version = "latest",
  cache_dir = interaction_database_cache_dir(),
  build_if_missing = FALSE,
  ...
) {
  database <- match.arg(database)
  assert_single_value(version, type = "string")
  assert_single_value(cache_dir, type = "string")
  assert_single_value(build_if_missing, type = "bool")

  path <- .interaction_database_rds_path(database, version, cache_dir)
  if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
    if (!isTRUE(build_if_missing)) {
      build_fun <- paste0("build_", database, "_database")
      cli::cli_abort(c(
        "x" = paste(
          "Interaction database {.val {database}} (version {.val {version}})",
          "not found at {.path {path}}."
        ),
        "i" = paste(
          "Run {.code load_interaction_database(\"{database}\", build_if_missing = TRUE)}",
          "or {.code {build_fun}()}."
        )
      ))
    }

    cli::cli_inform(c(
      "i" = paste(
        "Building interaction database {.val {database}}",
        "(version {.val {version}}) under {.path {cache_dir}}."
      )
    ))
    builder_args <- list(cache_dir = cache_dir, ...)
    if (!identical(version, "latest") && !"version" %in% names(builder_args)) {
      builder_args$version <- version
    }
    switch(
      database,
      string = do.call(build_string_database, builder_args),
      biogrid = do.call(build_biogrid_database, builder_args),
      corum = do.call(build_corum_database, builder_args),
      omnipath = do.call(build_omnipath_database, builder_args),
      alphafold = do.call(build_alphafold_database, builder_args)
    )
    if (!file.exists(path) || !isTRUE(file.info(path)$size > 0)) {
      cli::cli_abort(c(
        "x" = paste(
          "Build finished but {.val {database}} (version {.val {version}})",
          "is still missing at {.path {path}}."
        ),
        "i" = paste(
          "The builder may have written a different version label;",
          "check {.path {cache_dir}} or pass the concrete {.arg version}."
        )
      ))
    }
  }

  obj <- readRDS(path)
  assert_class(obj, "list")
  if (is.null(obj$edges)) {
    cli::cli_abort("Invalid interaction database object at {.path {path}}.")
  }
  assert_class(obj$edges, c("data.frame", "tbl_df"))
  required_cols <- c("uniprot_a", "uniprot_b")
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
  additional_columns <- obj$meta$additional_columns %||% character()
  if (!is.character(additional_columns)) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database meta at {.path {path}}.",
      "i" = "{.field additional_columns} must be a character vector."
    ))
  }
  missing_score_cols <- setdiff(score_columns, names(obj$edges))
  if (length(missing_score_cols) > 0) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database edges at {.path {path}}.",
      "i" = "Missing declared score columns: {.val {missing_score_cols}}"
    ))
  }
  missing_add_cols <- setdiff(additional_columns, names(obj$edges))
  if (length(missing_add_cols) > 0) {
    cli::cli_abort(c(
      "x" = "Invalid interaction database edges at {.path {path}}.",
      "i" = "Missing declared additional columns: {.val {missing_add_cols}}"
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
  obj$meta$additional_columns <- additional_columns
  obj$meta$resource <- obj$meta$resource %||% obj$meta$database %||% database
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
#' Maps panel markers to UniProt accessions, loads an interaction database,
#' and returns undirected marker pairs with a registered database entry after
#' optional score / network filters.
#' UniProt self-interactions are returned only as marker self-pairs. When
#' multiple marker names map to the same UniProt accession, a self-interaction
#' does not create interactions between those different marker names.
#' Conversely, interactions between different UniProt accessions are not
#' returned as marker self-pairs. Marker names are ordered alphabetically, and
#' the corresponding UniProt accessions are reordered with them.
#' When multiple UniProt edges collapse to the same marker pair, the row with
#' the highest score envelope is kept so UniProt IDs and additional columns
#' stay aligned with the reported scores. For databases without score columns,
#' the first remaining edge for that marker pair is kept.
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
#' @param string_network For STRING: \code{"physical"} or \code{"full"}
#'   (filters the \code{network} additional column).
#' @param cache_dir Interaction database cache directory.
#' @param version Database version label (\code{"latest"} or a built version).
#' @param build_if_missing Passed to \code{\link{load_interaction_database}}.
#'   Defaults to \code{FALSE}.
#' @return A tibble of panel edges with \code{marker_1}, \code{marker_2},
#'   \code{uniprot_a}, \code{uniprot_b}, native score columns, and any
#'   \code{additional_columns} from the database. Pass
#'   \code{marker_1}/\code{marker_2} columns to
#'   \code{ColocalizationHeatmap(highlight_pairs = ...)}.
#'
#' @section Score filtering:
#' Edge RDS files keep native score columns rather than a single synthetic
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
#'   score_min = c(combined_score = 400),
#'   build_if_missing = TRUE
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
  version = "latest",
  build_if_missing = FALSE
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
  assert_single_value(build_if_missing, type = "bool")

  db <- load_interaction_database(
    database = database,
    version = version,
    cache_dir = cache_dir,
    build_if_missing = build_if_missing
  )
  edges <- db$edges
  score_cols <- db$meta$score_columns %||% character()
  add_cols <- db$meta$additional_columns %||% character()

  # Zero-row result matching the loaded database score / additional schema.
  empty_panel_interactions <- function() {
    empty_like <- function(cols) {
      lapply(cols, function(col) {
        if (col %in% names(edges)) {
          x <- edges[[col]]
          if (is.numeric(x)) {
            return(numeric())
          }
          if (is.logical(x)) {
            return(logical())
          }
        }
        character()
      }) %>%
        setNames(cols) %>%
        as_tibble()
    }
    return(bind_cols(
      tibble(
        marker_1 = character(),
        marker_2 = character(),
        uniprot_a = character(),
        uniprot_b = character()
      ),
      empty_like(score_cols),
      empty_like(add_cols)
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
    if (!"network" %in% names(edges)) {
      cli::cli_abort(c(
        "x" = "STRING edge table lacks a {.val network} column.",
        "i" = "Rebuild the STRING database with the current builders."
      ))
    }
    edges <- edges %>%
      filter(network == string_network)
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
  joined <- edges %>%
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
      all_of(add_cols)
    )

  # Collapse duplicate marker pairs to one representative edge. When scores
  # exist, keep the row with the highest score envelope so UniProt / extras
  # stay aligned with the reported scores; otherwise keep the first edge.
  if (length(score_cols) > 0 && nrow(joined) > 0) {
    score_mat <- as.matrix(joined[score_cols])
    joined$.score_rank <- apply(score_mat, 1, function(r) {
      if (all(is.na(r))) {
        -Inf
      } else {
        max(r, na.rm = TRUE)
      }
    })
    out <- joined %>%
      group_by(marker_1, marker_2) %>%
      slice_max(order_by = .score_rank, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(-.score_rank)
  } else {
    out <- joined %>%
      group_by(marker_1, marker_2) %>%
      slice_head(n = 1) %>%
      ungroup()
  }
  return(out)
}

#' Build STRING physical / full link database (human)
#'
#' Maintainer helper: downloads STRING release files into ephemeral staging
#' (or reads an existing \code{raw_dir}), writes an edge RDS via
#' \code{\link{save_interaction_database}}, then deletes owned staging on
#' success.
#'
#' @param raw_dir Directory for STRING download files (links + aliases).
#'   When \code{NULL} (default), files are downloaded into a temporary staging
#'   directory that is removed after a successful build. Caller-supplied paths
#'   are left untouched.
#' @param version STRING version label (default \code{"12.0"}).
#' @param cache_dir Output interaction database cache.
#' @param species NCBI taxon (default human 9606).
#' @param include_full If \code{TRUE}, also include the full STRING network.
#' @return Path to saved RDS.
#' @export
build_string_database <- function(
  raw_dir = NULL,
  version = "12.0",
  cache_dir = interaction_database_cache_dir(),
  species = 9606,
  include_full = FALSE
) {
  rlang::check_installed("data.table")
  owned_staging <- is.null(raw_dir)
  if (owned_staging) {
    raw_dir <- .interaction_database_staging_dir("string")
  } else {
    assert_single_value(raw_dir, type = "string")
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  }
  success <- FALSE
  on.exit(
    .finalize_interaction_database_staging(raw_dir, owned_staging, success),
    add = TRUE
  )

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
      rename(uniprot_b = uniprot) %>%
      mutate(network = network_label)
    normalise_interaction_edges(
      m2,
      a_col = "uniprot_a",
      b_col = "uniprot_b",
      score_cols = "combined_score",
      additional_cols = "network"
    )
  }

  phys <- read_links(physical_gz, "physical")
  full <- if (isTRUE(include_full)) read_links(full_gz, "full") else NULL
  if (is.null(phys) && is.null(full)) {
    cli::cli_abort(
      "No STRING link files found (expected physical and/or full links .txt.gz)."
    )
  }
  edges <- bind_rows(phys, full)

  path <- save_interaction_database(
    edges = edges,
    database = "string",
    version = version,
    cache_dir = cache_dir,
    score_columns = "combined_score",
    additional_columns = "network",
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
  )
  success <- TRUE
  return(path)
}

#' Build BioGRID physical interaction database (human)
#'
#' Downloads the BioGRID multi-validated physical (MV-Physical) tab3 dump for
#' a real release id. \code{version = "latest"} fetches
#' \code{BIOGRID-MV-Physical-LATEST.tab3.zip} and resolves the concrete release
#' from the extracted \code{BIOGRID-MV-Physical-X.Y.Z.tab3.txt} filename.
#' A concrete \code{version} such as \code{"5.0.259"} downloads that release
#' from the BioGRID Release Archive.
#'
#' The edge RDS is always saved under the resolved BioGRID release id (and
#' copied to \code{biogrid_latest.rds}). Owned staging downloads are deleted
#' after a successful build.
#'
#' @param raw_file Optional path to a local
#'   \code{BIOGRID-MV-Physical-X.Y.Z.tab3.txt} file. When \code{NULL}, the
#'   release is downloaded into ephemeral staging.
#' @param version BioGRID release id (\code{X.Y.Z}) or \code{"latest"}.
#' @param cache_dir Output interaction database cache.
#' @return Path to saved RDS.
#' @export
build_biogrid_database <- function(
  raw_file = NULL,
  version = "latest",
  cache_dir = interaction_database_cache_dir()
) {
  rlang::check_installed("data.table")
  assert_single_value(version, type = "string")
  if (!identical(version, "latest") &&
    !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version)) {
    cli::cli_abort(c(
      "x" = "{.arg version} must be {.val latest} or a BioGRID release id like {.val 5.0.259}.",
      "i" = "Got {.val {version}}."
    ))
  }

  owned_staging <- is.null(raw_file)
  raw_dir <- if (owned_staging) {
    .interaction_database_staging_dir("biogrid")
  } else {
    NULL
  }
  success <- FALSE
  on.exit(
    .finalize_interaction_database_staging(raw_dir, owned_staging, success),
    add = TRUE
  )
  tab3_name_re <- "^BIOGRID-MV-Physical-([0-9]+\\.[0-9]+\\.[0-9]+)\\.tab3\\.txt$"

  version_from_tab3_name <- function(path) {
    bn <- basename(path)
    if (!grepl(tab3_name_re, bn, ignore.case = TRUE)) {
      return(NA_character_)
    }
    return(sub(tab3_name_re, "\\1", bn, ignore.case = TRUE))
  }

  if (!is.null(raw_file)) {
    assert_single_value(raw_file, type = "string")
    if (!file.exists(raw_file)) {
      cli::cli_abort("BioGRID raw file not found: {.path {raw_file}}.")
    }
    parsed <- version_from_tab3_name(raw_file)
    if (identical(version, "latest")) {
      if (is.na(parsed)) {
        cli::cli_abort(c(
          "x" = "When {.arg version} is {.val latest}, {.arg raw_file} must be named like
          {.val BIOGRID-MV-Physical-5.0.259.tab3.txt}.",
          "i" = "Got {.path {raw_file}}."
        ))
      }
      version <- parsed
    } else if (!is.na(parsed) && !identical(parsed, version)) {
      cli::cli_abort(c(
        "x" = "{.arg raw_file} release {.val {parsed}} does not match {.arg version} {.val {version}}.",
        "i" = "Path: {.path {raw_file}}"
      ))
    }
  } else {
    if (identical(version, "latest")) {
      zip_path <- file.path(raw_dir, "BIOGRID-MV-Physical-LATEST.tab3.zip")
      .download_if_missing(
        paste0(
          "https://downloads.thebiogrid.org/Download/BioGRID/",
          "Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
        ),
        zip_path
      )
      utils::unzip(zip_path, exdir = raw_dir)
      # Keep zip member paths (not basename) so nested archive layouts resolve.
      members <- utils::unzip(zip_path, list = TRUE)$Name
      hit <- members[grepl(tab3_name_re, basename(members), ignore.case = TRUE)]
      if (length(hit) != 1) {
        cli::cli_abort(c(
          "x" = "Expected exactly one versioned MV-Physical tab3 in {.path {zip_path}}.",
          "i" = "Found {length(hit)} matching member{?s}."
        ))
      }
      raw_file <- file.path(raw_dir, hit[[1]])
      version <- version_from_tab3_name(raw_file)
      if (is.na(version)) {
        cli::cli_abort(
          "Could not parse BioGRID release id from {.path {raw_file}}."
        )
      }
    } else {
      raw_file <- file.path(
        raw_dir,
        paste0("BIOGRID-MV-Physical-", version, ".tab3.txt")
      )
      if (!file.exists(raw_file) || !isTRUE(file.info(raw_file)$size > 0)) {
        zip_path <- file.path(
          raw_dir,
          paste0("BIOGRID-MV-Physical-", version, ".tab3.zip")
        )
        .download_if_missing(
          paste0(
            "https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/",
            "BIOGRID-", version, "/BIOGRID-MV-Physical-", version, ".tab3.zip"
          ),
          zip_path
        )
        utils::unzip(zip_path, exdir = raw_dir)
      }
    }
    if (!file.exists(raw_file) || !isTRUE(file.info(raw_file)$size > 0)) {
      cli::cli_abort("Missing BioGRID MV-Physical file {.path {raw_file}}.")
    }
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
  evidence <- if (!is.na(exp_col)) as.character(dt[[exp_col]]) else NULL
  pairs <- data.frame(
    uniprot_a = first_ac(dt[[a_col]]),
    uniprot_b = first_ac(dt[[b_col]]),
    stringsAsFactors = FALSE
  )
  additional_cols <- character()
  if (!is.null(evidence) && !is.na(exp_col)) {
    pairs[[exp_col]] <- evidence
    additional_cols <- exp_col
  }
  pairs <- pairs[!is.na(pairs$uniprot_a) & !is.na(pairs$uniprot_b), , drop = FALSE]
  edges <- normalise_interaction_edges(
    pairs,
    a_col = "uniprot_a",
    b_col = "uniprot_b",
    additional_cols = if (length(additional_cols)) additional_cols else NULL
  )

  path <- save_interaction_database(
    edges = edges,
    database = "biogrid",
    version = version,
    cache_dir = cache_dir,
    score_columns = NULL,
    additional_columns = if (length(additional_cols)) additional_cols else NULL,
    source_url = paste0(
      "https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/",
      "BIOGRID-", version, "/"
    ),
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
  )
  success <- TRUE
  return(path)
}

#' Build CORUM co-membership database (human)
#'
#' Downloads OmniPath-served CORUM complexes when \code{corum_file} is missing.
#' Prefer an official CORUM file when present. Complexes with a single UniProt
#' accession become homomer edges (\code{uniprot_a == uniprot_b}). Owned staging
#' downloads are deleted after a successful build.
#'
#' @param corum_file Path to OmniPath CORUM complexes TSV or official
#'   coreComplexes. When \code{NULL} (default), the file is downloaded into
#'   ephemeral staging.
#' @param version Version label.
#' @param cache_dir Output interaction database cache.
#' @return Path to saved RDS.
#' @export
build_corum_database <- function(
  corum_file = NULL,
  version = "omnipath_corum",
  cache_dir = interaction_database_cache_dir()
) {
  rlang::check_installed("data.table")
  owned_staging <- is.null(corum_file)
  staging_dir <- NULL
  if (owned_staging) {
    staging_dir <- .interaction_database_staging_dir("corum")
    corum_file <- file.path(staging_dir, "omnipath_complexes_corum.tsv")
  } else {
    assert_single_value(corum_file, type = "string")
  }
  success <- FALSE
  on.exit(
    .finalize_interaction_database_staging(staging_dir, owned_staging, success),
    add = TRUE
  )
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

  pieces <- lapply(seq_len(nrow(dt)), function(i) {
    ids <- unique(trimws(unlist(strsplit(
      as.character(dt[[comp_col]][[i]]),
      "[_;,]"
    ))))
    ids <- ids[ids != "" & !is.na(ids)]
    if (length(ids) < 1) {
      return(NULL)
    }
    # Single-accession complexes (homomers) become U-U edges
    if (length(ids) == 1) {
      row <- tibble(uniprot_a = ids, uniprot_b = ids)
    } else {
      grid <- utils::combn(ids, 2)
      row <- tibble(
        uniprot_a = grid[1, ],
        uniprot_b = grid[2, ]
      )
    }
    if (!is.na(name_col)) {
      row[[name_col]] <- as.character(dt[[name_col]][[i]])
    }
    row
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) < 1) {
    pairs <- tibble(uniprot_a = character(), uniprot_b = character())
  } else {
    pairs <- bind_rows(pieces)
  }
  additional_cols <- if (!is.na(name_col) && name_col %in% names(pairs)) {
    name_col
  } else {
    NULL
  }
  edges <- normalise_interaction_edges(
    pairs,
    additional_cols = additional_cols
  )

  path <- save_interaction_database(
    edges = edges,
    database = "corum",
    version = version,
    cache_dir = cache_dir,
    score_columns = NULL,
    additional_columns = additional_cols,
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
  )
  success <- TRUE
  return(path)
}

#' Build OmniPath interaction database
#'
#' Downloads OmniPath interactions into ephemeral staging when
#' \code{interactions_file} is \code{NULL}. Owned staging is deleted after a
#' successful build.
#'
#' @param interactions_file Path to OmniPath interactions TSV. When
#'   \code{NULL}, downloaded into ephemeral staging (not for
#'   \code{license = "unknown"}, which requires a local file).
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
  owned_staging <- is.null(interactions_file)
  staging_dir <- NULL
  if (owned_staging) {
    if (identical(license, "unknown")) {
      cli::cli_abort(c(
        "x" = "Missing OmniPath interactions file.",
        "i" = paste(
          "Pass {.arg interactions_file} when {.arg license} is {.val unknown}."
        )
      ))
    }
    staging_dir <- .interaction_database_staging_dir("omnipath")
    interactions_file <- file.path(
      staging_dir,
      paste0("interactions_omnipath_", license, ".tsv")
    )
  } else {
    assert_single_value(interactions_file, type = "string")
  }
  success <- FALSE
  on.exit(
    .finalize_interaction_database_staging(staging_dir, owned_staging, success),
    add = TRUE
  )
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
  pairs <- data.frame(
    uniprot_a = as.character(dt[[src]]),
    uniprot_b = as.character(dt[[tgt]]),
    stringsAsFactors = FALSE
  )
  additional_cols <- intersect(c("sources", "references"), colnames(dt))
  for (col in additional_cols) {
    pairs[[col]] <- as.character(dt[[col]])
  }
  edges <- normalise_interaction_edges(
    pairs,
    additional_cols = if (length(additional_cols)) additional_cols else NULL
  )

  path <- save_interaction_database(
    edges = edges,
    database = "omnipath",
    version = paste0(version, "_", license),
    cache_dir = cache_dir,
    score_columns = NULL,
    additional_columns = if (length(additional_cols)) additional_cols else NULL,
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
      "resources identified in the sources/references fields"
    )
  )
  success <- TRUE
  return(path)
}

#' Build AlphaFold DB complex database
#'
#' Combines heterodimer and homodimer predictions. When no local CSVs are
#' present, downloads both NVIDIA/AFDB metadata tables from EBI FTP into
#' ephemeral staging (heterodimer ~2 GB, homodimer ~6 GB), writes an edge RDS
#' with all parsed human edges and native scores, then deletes owned staging
#' on success. Apply confidence cuts at query time with
#' \code{\link{extract_panel_interactions}} (\code{score_min} /
#' \code{score_max}).
#'
#' Homodimer metadata uses columns \code{uniprotAccession} and \code{taxId};
#' heterodimers use \code{uniprot_ac_1}/\code{uniprot_ac_2} and
#' \code{tax_id_1}/\code{tax_id_2}.
#'
#' @param heterodimer_file Optional local heterodimer metadata CSV. Ignored when
#'   missing; panel CSVs under \code{raw_dir} are tried next.
#' @param raw_dir Directory for AFDB panel CSVs and/or official dumps. When
#'   \code{NULL} (default), a temporary staging directory is used and removed
#'   after a successful build. Caller-supplied paths are left untouched.
#' @param version Version label.
#' @param cache_dir Output interaction database cache.
#' @return Path to saved RDS.
#' @export
build_alphafold_database <- function(
  heterodimer_file = NULL,
  raw_dir = NULL,
  version = "afdb_nvda",
  cache_dir = interaction_database_cache_dir()
) {
  owned_staging <- is.null(raw_dir)
  if (owned_staging) {
    raw_dir <- .interaction_database_staging_dir("alphafold")
  } else {
    assert_single_value(raw_dir, type = "string")
    dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  }
  success <- FALSE
  on.exit(
    .finalize_interaction_database_staging(raw_dir, owned_staging, success),
    add = TRUE
  )

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

  panel_files <- c(
    if (!is.null(heterodimer_file)) heterodimer_file else character(),
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    file.path(raw_dir, "afdb_heterodimers_panel_pairs.csv"),
    file.path(raw_dir, "afdb_homodimers_human_panel.csv")
  )
  panel_files <- unique(panel_files[file.exists(panel_files)])
  has_panel <- length(panel_files) > 0
  has_both_official <- file.exists(hetero_dest) &&
    isTRUE(file.info(hetero_dest)$size > 0) &&
    file.exists(homo_dest) &&
    isTRUE(file.info(homo_dest)$size > 0)

  # Panel/fixture CSVs: use them and do not pull multi-GB FTP dumps.
  # Otherwise ensure both official NVIDIA dumps are present (a single leftover
  # dump must not skip its counterpart after a partial download).
  if (!has_panel && !has_both_official) {
    cli::cli_inform(c(
      "i" = paste(
        "Downloading NVIDIA AFDB heterodimer (~2 GB) and/or",
        "homodimer (~6 GB) metadata as needed."
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
  }

  official <- c(hetero_dest, homo_dest)
  official_ok <- vapply(official, function(path) {
    file.exists(path) && isTRUE(file.info(path)$size > 0)
  }, logical(1))
  files <- unique(c(panel_files, official[official_ok]))

  pieces <- lapply(files, parse_afdb_file)
  pairs <- bind_rows(pieces)
  if (nrow(pairs) < 1) {
    cli::cli_abort("No usable AFDB pairs found in {.path {raw_dir}}.")
  }
  edges <- normalise_interaction_edges(
    pairs,
    score_cols = c("ipSAE", "pDockQ2")
  )

  path <- save_interaction_database(
    edges = edges,
    database = "alphafold",
    version = version,
    cache_dir = cache_dir,
    score_columns = c("ipSAE", "pDockQ2"),
    additional_columns = NULL,
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
  )
  success <- TRUE
  return(path)
}

#' Build all five interaction databases
#'
#' Maintainer helper that runs each \code{build_*_database()} writer. Missing
#' raw dumps are downloaded into ephemeral staging and removed after each
#' successful build. STRING is built with both physical and full networks
#' (\code{include_full = TRUE}). OmniPath defaults to the commercial license
#' filter.
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
