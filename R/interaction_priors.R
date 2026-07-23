#' Interaction prior helpers (slim UniProt edge tables)
#'
#' Canonical edge schema:
#' \code{uniprot_a}, \code{uniprot_b} (undirected, a <= b), optional
#' \code{score}, \code{evidence}, plus \code{resource}, \code{resource_version},
#' \code{built_at}.

.interaction_prior_databases <- c(
  "string", "biogrid", "corum", "omnipath", "alphafold"
)

#' Default cache directory for interaction priors
#'
#' @param cache_dir Optional override. If \code{NULL}, uses
#'   \code{results/db_cache/priors} under the project root when \code{here} is
#'   available, otherwise \code{tools::R_user_dir("pixelatorR", "cache")/priors}.
#' @return Character path.
#' @export
InteractionPriorCacheDir <- function(cache_dir = NULL) {
  if (!is.null(cache_dir)) {
    return(cache_dir)
  }
  if (requireNamespace("here", quietly = TRUE)) {
    return(here::here("results", "db_cache", "priors"))
  }
  file.path(tools::R_user_dir("pixelatorR", which = "cache"), "priors")
}

#' Make an undirected UniProt pair key
#'
#' @param a,b Character vectors of UniProt accessions.
#' @return Character vector of keys \code{"A|B"} with A <= B.
#' @export
UndirectedUniprotPairKey <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  paste(pmin(a, b), pmax(a, b), sep = "|")
}

#' Normalise an edge table to the canonical prior schema
#'
#' @param edges Data frame with at least two UniProt columns.
#' @param a_col,b_col Column names for the two ends.
#' @param score_col Optional score column.
#' @param evidence_col Optional evidence column.
#' @param resource Resource name.
#' @param resource_version Version string.
#' @return Tibble with canonical columns.
#' @export
NormaliseInteractionPrior <- function(
  edges,
  a_col = "uniprot_a",
  b_col = "uniprot_b",
  score_col = NULL,
  evidence_col = NULL,
  resource,
  resource_version
) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required", call. = FALSE)
  }
  edges <- as.data.frame(edges)
  if (!all(c(a_col, b_col) %in% names(edges))) {
    stop("Missing UniProt columns '", a_col, "' / '", b_col, "'", call. = FALSE)
  }
  a <- trimws(as.character(edges[[a_col]]))
  b <- trimws(as.character(edges[[b_col]]))
  keep <- !is.na(a) & !is.na(b) & a != "" & b != ""
  a <- a[keep]
  b <- b[keep]
  out <- tibble::tibble(
    uniprot_a = pmin(a, b),
    uniprot_b = pmax(a, b),
    score = if (!is.null(score_col) && score_col %in% names(edges)) {
      as.numeric(edges[[score_col]][keep])
    } else {
      NA_real_
    },
    evidence = if (!is.null(evidence_col) && evidence_col %in% names(edges)) {
      as.character(edges[[evidence_col]][keep])
    } else {
      NA_character_
    },
    resource = as.character(resource),
    resource_version = as.character(resource_version),
    built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
  out <- unique(out)
  out
}

.prior_rds_path <- function(database, version, cache_dir) {
  file.path(cache_dir, paste0(database, "_", version, ".rds"))
}

.prior_manifest_path <- function(cache_dir) {
  file.path(cache_dir, "MANIFEST.json")
}

.update_prior_manifest <- function(cache_dir, entry) {
  path <- .prior_manifest_path(cache_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  man <- list()
  if (file.exists(path)) {
    man <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
  man[[entry$database]] <- entry
  jsonlite::write_json(man, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}

#' Write a slim interaction prior RDS + manifest entry
#'
#' @param edges Canonical edge tibble.
#' @param database Database key.
#' @param version Version label used in the filename.
#' @param cache_dir Cache directory.
#' @param source_url Optional provenance URL.
#' @param license Optional license string.
#' @param citation Optional citation string.
#' @return Path to the written RDS.
#' @export
SaveInteractionPrior <- function(
  edges,
  database,
  version,
  cache_dir = InteractionPriorCacheDir(),
  source_url = NULL,
  license = NULL,
  citation = NULL
) {
  database <- match.arg(database, .interaction_prior_databases)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- .prior_rds_path(database, version, cache_dir)
  meta <- list(
    database = database,
    version = version,
    n_edges = nrow(edges),
    source_url = source_url,
    license = license,
    citation = citation,
    built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    file = basename(path)
  )
  saveRDS(list(edges = edges, meta = meta), path)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    .update_prior_manifest(cache_dir, meta)
  }
  # also write/update "latest" pointer copy
  latest <- .prior_rds_path(database, "latest", cache_dir)
  file.copy(path, latest, overwrite = TRUE)
  path
}

#' Load a slim interaction prior
#'
#' @param database One of string, biogrid, corum, omnipath, alphafold.
#' @param version Version label or \code{"latest"}.
#' @param cache_dir Cache directory.
#' @return List with \code{edges} (tibble) and \code{meta}.
#' @export
LoadInteractionPrior <- function(
  database = c("string", "biogrid", "corum", "omnipath", "alphafold"),
  version = "latest",
  cache_dir = InteractionPriorCacheDir()
) {
  database <- match.arg(database)
  path <- .prior_rds_path(database, version, cache_dir)
  if (!file.exists(path)) {
    stop(
      "Interaction prior not found: ", path, "\n",
      "Build it with Build",
      paste0(toupper(substring(database, 1, 1)), substring(database, 2)),
      "Prior() (see R/BuildInteractionPriors.R).",
      call. = FALSE
    )
  }
  obj <- readRDS(path)
  if (!is.list(obj) || is.null(obj$edges)) {
    stop("Invalid prior object at ", path, call. = FALSE)
  }
  obj
}
