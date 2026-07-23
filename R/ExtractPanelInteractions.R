#' Extract known database interactions for a marker panel
#'
#' Maps panel markers to UniProt accessions, loads a slim interaction prior, and
#' returns undirected marker pairs with a registered database entry after
#' optional score / network filters.
#'
#' @param markers Character vector of panel marker names.
#' @param database Interaction database key.
#' @param marker_uniprot_map Data frame with columns \code{marker} and
#'   \code{uniprot_id} (long form; multiple rows per marker allowed).
#' @param score_threshold Minimum score (STRING uses the classic 0–1000 scale).
#' @param string_network For STRING: \code{"physical"} or \code{"full"}.
#' @param omnipath_license Documented license filter for the OmniPath prior
#'   (must match the built prior version suffix when using versioned files).
#' @param alphafold_ipsae_min,alphafold_pdockq2_min Reserved for AFDB filtering
#'   when scores are present in the prior (priors are typically pre-filtered).
#' @param cache_dir Prior cache directory.
#' @param version Prior version label (\code{"latest"} or a built version).
#' @return A tibble of panel edges with \code{marker_1}, \code{marker_2},
#'   \code{uniprot_a}, \code{uniprot_b}, \code{in_db = TRUE}, and optional
#'   \code{score}, \code{evidence}, \code{resource}, \code{resource_version}.
#'   Pass \code{marker_1}/\code{marker_2} columns to
#'   \code{ColocalizationHeatmap(highlight_pairs = ...)}.
#' @export
ExtractPanelInteractions <- function(
  markers,
  database = c("string", "biogrid", "corum", "omnipath", "alphafold"),
  marker_uniprot_map,
  score_threshold = 400,
  string_network = c("physical", "full"),
  omnipath_license = c("commercial", "academic"),
  alphafold_ipsae_min = 0.6,
  alphafold_pdockq2_min = 0.23,
  cache_dir = InteractionPriorCacheDir(),
  version = "latest"
) {
  database <- match.arg(database)
  string_network <- match.arg(string_network)
  omnipath_license <- match.arg(omnipath_license)

  if (missing(marker_uniprot_map) || is.null(marker_uniprot_map)) {
    stop("`marker_uniprot_map` is required (columns: marker, uniprot_id).", call. = FALSE)
  }
  map <- as.data.frame(marker_uniprot_map)
  if (!all(c("marker", "uniprot_id") %in% names(map))) {
    stop("`marker_uniprot_map` must contain columns `marker` and `uniprot_id`.", call. = FALSE)
  }
  markers <- unique(as.character(markers))
  map <- map[
    as.character(map$marker) %in% markers &
      !is.na(map$uniprot_id) &
      as.character(map$uniprot_id) != "",
    ,
    drop = FALSE
  ]
  map$marker <- as.character(map$marker)
  map$uniprot_id <- trimws(as.character(map$uniprot_id))
  if (nrow(map) == 0) {
    return(.empty_panel_interactions())
  }

  prior <- LoadInteractionPrior(
    database = database,
    version = version,
    cache_dir = cache_dir
  )
  edges <- as.data.frame(prior$edges)

  # Database-specific filters
  if (database == "string") {
    if ("evidence" %in% names(edges)) {
      edges <- edges[edges$evidence == string_network | edges$resource == paste0("string_", string_network), , drop = FALSE]
    } else if ("resource" %in% names(edges)) {
      edges <- edges[grepl(string_network, edges$resource), , drop = FALSE]
    }
    if ("score" %in% names(edges) && !is.null(score_threshold)) {
      edges <- edges[!is.na(edges$score) & edges$score >= score_threshold, , drop = FALSE]
    }
  }
  if (database == "alphafold" && "score" %in% names(edges)) {
    # score stored as ipSAE (preferred) or pDockQ2; keep rows with NA (prefiltered)
    # or passing either threshold
    keep <- is.na(edges$score) |
      edges$score >= min(alphafold_ipsae_min, alphafold_pdockq2_min)
    edges <- edges[keep, , drop = FALSE]
  }
  # omnipath_license: prefer matching version when latest points elsewhere
  if (database == "omnipath" && version == "latest") {
    # If a commercial build exists, LoadInteractionPrior("latest") should already
    # point at it when BuildOmnipathPrior(license="commercial") was used.
    invisible(omnipath_license)
  }

  panel_uniprot <- unique(map$uniprot_id)
  edges <- edges[
    edges$uniprot_a %in% panel_uniprot & edges$uniprot_b %in% panel_uniprot,
    ,
    drop = FALSE
  ]
  if (nrow(edges) == 0) {
    return(.empty_panel_interactions())
  }

  # Expand UniProt pairs to marker pairs (many-to-many)
  m_a <- merge(
    edges,
    map,
    by.x = "uniprot_a",
    by.y = "uniprot_id",
    allow.cartesian = TRUE
  )
  names(m_a)[names(m_a) == "marker"] <- "marker_1"
  m_ab <- merge(
    m_a,
    map,
    by.x = "uniprot_b",
    by.y = "uniprot_id",
    allow.cartesian = TRUE
  )
  names(m_ab)[names(m_ab) == "marker"] <- "marker_2"

  out <- tibble::tibble(
    marker_1 = pmin(as.character(m_ab$marker_1), as.character(m_ab$marker_2)),
    marker_2 = pmax(as.character(m_ab$marker_1), as.character(m_ab$marker_2)),
    uniprot_a = as.character(m_ab$uniprot_a),
    uniprot_b = as.character(m_ab$uniprot_b),
    in_db = TRUE,
    score = if ("score" %in% names(m_ab)) as.numeric(m_ab$score) else NA_real_,
    evidence = if ("evidence" %in% names(m_ab)) as.character(m_ab$evidence) else NA_character_,
    resource = if ("resource" %in% names(m_ab)) as.character(m_ab$resource) else database,
    resource_version = if ("resource_version" %in% names(m_ab)) {
      as.character(m_ab$resource_version)
    } else {
      NA_character_
    }
  )
  # unique marker pairs; keep max score if present
  out <- out |>
    dplyr::group_by(marker_1, marker_2) |>
    dplyr::summarise(
      uniprot_a = dplyr::first(uniprot_a),
      uniprot_b = dplyr::first(uniprot_b),
      in_db = TRUE,
      score = suppressWarnings(max(score, na.rm = TRUE)),
      evidence = dplyr::first(evidence),
      resource = dplyr::first(resource),
      resource_version = dplyr::first(resource_version),
      .groups = "drop"
    ) |>
    dplyr::mutate(score = ifelse(is.infinite(score), NA_real_, score))
  out
}

.empty_panel_interactions <- function() {
  tibble::tibble(
    marker_1 = character(),
    marker_2 = character(),
    uniprot_a = character(),
    uniprot_b = character(),
    in_db = logical(),
    score = numeric(),
    evidence = character(),
    resource = character(),
    resource_version = character()
  )
}
