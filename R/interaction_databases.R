# Canonical edge schema: uniprot_a, uniprot_b (undirected, a <= b), optional
# score, evidence, plus resource, resource_version, built_at.

.interaction_databases <- c(
  "string", "biogrid", "corum", "omnipath", "alphafold"
)

.default_raw_cache <- function() {
  file.path(tools::R_user_dir("pixelatorR", which = "cache"), "db_raw")
}

.interaction_database_rds_path <- function(database, version, cache_dir) {
  file.path(cache_dir, paste0(database, "_", version, ".rds"))
}

#' Download a file if it is missing or empty
#' @noRd
.download_if_missing <- function(url, dest_file) {
  if (file.exists(dest_file) && isTRUE(file.info(dest_file)$size > 0)) {
    return(dest_file)
  }
  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
  cli::cli_inform(c("i" = "Downloading {.url {url}}"))
  status <- tryCatch(
    utils::download.file(url, destfile = dest_file, mode = "wb", quiet = TRUE),
    error = function(e) {
      cli::cli_abort(c(
        "x" = "Failed to download {.url {url}}",
        "i" = conditionMessage(e)
      ))
    }
  )
  if (!identical(as.integer(status), 0L) ||
      !file.exists(dest_file) ||
      !isTRUE(file.info(dest_file)$size > 0)) {
    cli::cli_abort("Download failed or produced an empty file: {.url {url}}")
  }
  return(dest_file)
}

#' Download a zip archive and return the first extracted file matching pattern
#' @noRd
.download_zip_extract <- function(url, dest_dir, pattern) {
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

.string_download_url <- function(fname, version, network = c("aliases", "physical", "full")) {
  network <- match.arg(network)
  prefix <- switch(
    network,
    aliases = paste0("protein.aliases.v", version),
    physical = paste0("protein.physical.links.v", version),
    full = paste0("protein.links.v", version)
  )
  paste0("https://stringdb-downloads.org/download/", prefix, "/", fname)
}

.ensure_string_file <- function(raw_dir, fname, version, network) {
  path <- file.path(raw_dir, fname)
  if (file.exists(path) && isTRUE(file.info(path)$size > 0)) {
    return(path)
  }
  return(.download_if_missing(
    .string_download_url(fname, version, network),
    path
  ))
}

.empty_panel_interactions <- function() {
  tibble(
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
#' @param score_col Optional score column.
#' @param evidence_col Optional evidence column.
#' @param resource Resource name.
#' @param resource_version Version string.
#' @return Tibble with canonical columns.
#' @export
normalise_interaction_edges <- function(
  edges,
  a_col = "uniprot_a",
  b_col = "uniprot_b",
  score_col = NULL,
  evidence_col = NULL,
  resource,
  resource_version
) {
  assert_class(edges, c("data.frame", "tbl_df"))
  assert_single_value(a_col, type = "string")
  assert_single_value(b_col, type = "string")
  assert_single_value(score_col, type = "string", allow_null = TRUE)
  assert_single_value(evidence_col, type = "string", allow_null = TRUE)
  assert_single_value(resource, type = "string")
  assert_single_value(resource_version, type = "string")
  assert_col_in_data(a_col, edges)
  assert_col_in_data(b_col, edges)

  a <- trimws(as.character(edges[[a_col]]))
  b <- trimws(as.character(edges[[b_col]]))
  keep <- !is.na(a) & !is.na(b) & a != "" & b != ""
  a <- a[keep]
  b <- b[keep]
  return(unique(tibble(
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
    resource = resource,
    resource_version = resource_version,
    built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )))
}

#' Write a slim interaction database RDS
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
save_interaction_database <- function(
  edges,
  database = .interaction_databases,
  version,
  cache_dir = interaction_database_cache_dir(),
  source_url = NULL,
  license = NULL,
  citation = NULL
) {
  database <- match.arg(database, choices = .interaction_databases)
  assert_class(edges, c("data.frame", "tbl_df"))
  assert_single_value(version, type = "string")
  assert_single_value(cache_dir, type = "string")
  assert_single_value(source_url, type = "string", allow_null = TRUE)
  assert_single_value(license, type = "string", allow_null = TRUE)
  assert_single_value(citation, type = "string", allow_null = TRUE)

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- .interaction_database_rds_path(database, version, cache_dir)
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
  latest <- .interaction_database_rds_path(database, "latest", cache_dir)
  file.copy(path, latest, overwrite = TRUE)
  return(path)
}

#' Load a slim interaction database
#'
#' @param database One of \code{.interaction_databases}.
#' @param version Version label or \code{"latest"}.
#' @param cache_dir Cache directory.
#' @return List with \code{edges} (tibble) and \code{meta}.
#' @export
load_interaction_database <- function(
  database = .interaction_databases,
  version = "latest",
  cache_dir = interaction_database_cache_dir()
) {
  database <- match.arg(database, choices = .interaction_databases)
  assert_single_value(version, type = "string")
  assert_single_value(cache_dir, type = "string")

  path <- .interaction_database_rds_path(database, version, cache_dir)
  assert_file_exists(path)
  obj <- readRDS(path)
  assert_class(obj, "list")
  if (is.null(obj$edges)) {
    cli::cli_abort("Invalid interaction database object at {.path {path}}.")
  }
  return(obj)
}

#' Extract known database interactions for a marker panel
#'
#' Maps panel markers to UniProt accessions, loads a slim interaction database,
#' and returns undirected marker pairs with a registered database entry after
#' optional score / network filters.
#'
#' @param markers Character vector of panel marker names.
#' @param database Interaction database key.
#' @param marker_uniprot_map Data frame/tibble with columns \code{marker} and
#'   \code{uniprot_id} (long form; multiple rows per marker allowed).
#' @param score_threshold Minimum score (STRING uses the classic 0–1000 scale).
#' @param string_network For STRING: \code{"physical"} or \code{"full"}.
#' @param cache_dir Interaction database cache directory.
#' @param version Database version label (\code{"latest"} or a built version).
#' @return A tibble of panel edges with \code{marker_1}, \code{marker_2},
#'   \code{uniprot_a}, \code{uniprot_b}, \code{in_db = TRUE}, and optional
#'   \code{score}, \code{evidence}, \code{resource}, \code{resource_version}.
#'   Pass \code{marker_1}/\code{marker_2} columns to
#'   \code{ColocalizationHeatmap(highlight_pairs = ...)}.
#'
#' @examples
#' \dontrun{
#' highlight_pairs <- extract_panel_interactions(
#'   markers = c("CD3E", "CD4", "CD8A"),
#'   database = "string",
#'   marker_uniprot_map = marker_map,
#'   score_threshold = 400
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
  database = .interaction_databases,
  marker_uniprot_map,
  score_threshold = 400,
  string_network = c("physical", "full"),
  cache_dir = interaction_database_cache_dir(),
  version = "latest"
) {
  database <- match.arg(database, choices = .interaction_databases)
  string_network <- match.arg(string_network, choices = c("physical", "full"))
  assert_vector(markers, type = "character", n = 1)
  assert_class(marker_uniprot_map, c("data.frame", "tbl_df"))
  assert_col_in_data("marker", marker_uniprot_map)
  assert_col_in_data("uniprot_id", marker_uniprot_map)
  assert_col_class("marker", marker_uniprot_map, classes = "character")
  assert_col_class("uniprot_id", marker_uniprot_map, classes = "character")
  assert_single_value(score_threshold, type = "numeric", allow_null = TRUE)
  assert_single_value(cache_dir, type = "string")
  assert_single_value(version, type = "string")

  markers <- unique(markers)
  map <- marker_uniprot_map %>%
    filter(
      marker %in% markers,
      !is.na(uniprot_id),
      uniprot_id != ""
    ) %>%
    mutate(uniprot_id = trimws(uniprot_id))
  if (nrow(map) == 0) {
    return(.empty_panel_interactions())
  }

  db <- load_interaction_database(
    database = database,
    version = version,
    cache_dir = cache_dir
  )
  edges <- db$edges

  if (database == "string") {
    if ("evidence" %in% names(edges)) {
      edges <- edges[edges$evidence == string_network, , drop = FALSE]
    }
    if ("score" %in% names(edges) && !is.null(score_threshold)) {
      edges <- edges[
        !is.na(edges$score) & edges$score >= score_threshold,
        ,
        drop = FALSE
      ]
    }
  }

  panel_uniprot <- unique(map$uniprot_id)
  edges <- edges %>%
    filter(uniprot_a %in% panel_uniprot,
           uniprot_b %in% panel_uniprot)

  if (nrow(edges) == 0) {
    return(.empty_panel_interactions())
  }

  m_ab <- edges %>%
    left_join(map, by = c("uniprot_a" = "uniprot_id")) %>%
    rename(marker_1 = marker) %>%
    left_join(map, by = c("uniprot_b" = "uniprot_id")) %>%
    rename(marker_2 = marker)

  out <- tibble(
    marker_1 = pmin(m_ab$marker_1, m_ab$marker_2),
    marker_2 = pmax(m_ab$marker_1, m_ab$marker_2),
    uniprot_a = m_ab$uniprot_a,
    uniprot_b = m_ab$uniprot_b,
    in_db = TRUE,
    score = if ("score" %in% names(m_ab)) as.numeric(m_ab$score) else NA_real_,
    evidence = if ("evidence" %in% names(m_ab)) {
      as.character(m_ab$evidence)
    } else {
      NA_character_
    },
    resource = if ("resource" %in% names(m_ab)) {
      as.character(m_ab$resource)
    } else {
      database
    },
    resource_version = if ("resource_version" %in% names(m_ab)) {
      as.character(m_ab$resource_version)
    } else {
      NA_character_
    }
  ) %>%
    group_by(marker_1, marker_2) %>%
    summarise(
      uniprot_a = uniprot_a[1],
      uniprot_b = uniprot_b[1],
      in_db = TRUE,
      score = suppressWarnings(max(score, na.rm = TRUE)),
      evidence = evidence[1],
      resource = resource[1],
      resource_version = resource_version[1],
      .groups = "drop"
    ) %>%
    mutate(score = ifelse(is.infinite(score), NA_real_, score))
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

  aliases_fname <- sprintf("%s.protein.aliases.v%s.txt.gz", species, version)
  physical_fname <- sprintf(
    "%s.protein.physical.links.v%s.txt.gz", species, version
  )
  full_fname <- sprintf("%s.protein.links.v%s.txt.gz", species, version)

  aliases_gz <- .ensure_string_file(
    raw_dir, aliases_fname, version, "aliases"
  )
  physical_gz <- .ensure_string_file(
    raw_dir, physical_fname, version, "physical"
  )
  full_gz <- if (isTRUE(include_full)) {
    .ensure_string_file(raw_dir, full_fname, version, "full")
  } else {
    path <- file.path(raw_dir, full_fname)
    if (file.exists(path)) path else NA_character_
  }

  aliases <- data.table::fread(
    aliases_gz,
    header = TRUE,
    sep = "\t",
    quote = "",
    showProgress = FALSE
  )
  data.table::setnames(aliases, c("string_protein_id", "alias", "source"))
  .pick_uniprot <- function(aliases_vec) {
    swiss <- aliases_vec[grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", aliases_vec)]
    if (length(swiss)) {
      return(swiss[[1]])
    }
    aliases_vec[[1]]
  }
  up <- aliases[
    source == "UniProt_AC",
    .(uniprot = .pick_uniprot(alias)),
    by = string_protein_id
  ]

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
      left_join(as_tibble(up), by = c("protein1" = "string_protein_id")) %>%
      rename(uniprot_a = uniprot) %>%
      left_join(as_tibble(up), by = c("protein2" = "string_protein_id")) %>%
      rename(uniprot_b = uniprot)
    edges <- normalise_interaction_edges(
      m2,
      a_col = "uniprot_a",
      b_col = "uniprot_b",
      score_col = "combined_score",
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
    source_url = "https://string-db.org/cgi/download",
    license = "CC BY 4.0",
    citation = "Szklarczyk et al.; STRING database (see string-db.org)"
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
  if (is.null(raw_file)) {
    existing <- list.files(raw_dir, pattern = "[.]tab3[.]txt$", full.names = TRUE)
    raw_file <- if (length(existing) > 0) existing[[1]] else NA_character_
  }
  if (length(raw_file) != 1 || is.na(raw_file) || !file.exists(raw_file)) {
    raw_file <- .download_zip_extract(
      url = paste0(
        "https://downloads.thebiogrid.org/Download/BioGRID/",
        "Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"
      ),
      dest_dir = raw_dir,
      pattern = "[.]tab3[.]txt$"
    )
  }
  dt <- data.table::fread(raw_file, sep = "\t", quote = "", showProgress = FALSE)
  a_col <- grep(
    "SWISS-PROT.*Interactor A|SWISS-PROT Accessions Interactor A",
    colnames(dt),
    value = TRUE,
    ignore.case = TRUE
  )
  b_col <- grep(
    "SWISS-PROT.*Interactor B|SWISS-PROT Accessions Interactor B",
    colnames(dt),
    value = TRUE,
    ignore.case = TRUE
  )
  org_a <- grep("Organism ID Interactor A", colnames(dt), value = TRUE)[1]
  org_b <- grep("Organism ID Interactor B", colnames(dt), value = TRUE)[1]
  exp_col <- grep("^Experimental System$", colnames(dt), value = TRUE)[1]
  if (length(a_col) < 1 || length(b_col) < 1) {
    cli::cli_abort("Could not find SWISS-PROT columns in {.path {raw_file}}.")
  }
  a_col <- a_col[1]
  b_col <- b_col[1]

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
    source_url = "https://downloads.thebiogrid.org/",
    license = "MIT",
    citation = "Stark et al., Nucleic Acids Res. 2006; BioGRID"
  ))
}

#' Build CORUM co-membership database (human)
#'
#' Downloads OmniPath-served CORUM complexes when \code{corum_file} is missing.
#' Prefer an official CORUM file when present.
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
  comp_col <- intersect(
    c(
      "components",
      "subunits_uniprot_id",
      "subunits",
      "subunits(UniProt IDs)"
    ),
    colnames(dt)
  )
  if (length(comp_col) < 1) {
    cli::cli_abort("No components column in {.path {corum_file}}.")
  }
  comp_col <- comp_col[1]
  name_col <- intersect(c("name", "complex_name", "ComplexName"), colnames(dt))
  name_col <- if (length(name_col)) name_col[1] else NULL

  pairs_a <- character()
  pairs_b <- character()
  evid <- character()
  comps <- as.character(dt[[comp_col]])
  names_vec <- if (!is.null(name_col)) {
    as.character(dt[[name_col]])
  } else {
    rep(NA_character_, length(comps))
  }
  for (i in seq_along(comps)) {
    ids <- unique(trimws(unlist(strsplit(comps[i], "[_;,]"))))
    ids <- ids[ids != "" & !is.na(ids)]
    if (length(ids) < 2) {
      next
    }
    grid <- utils::combn(ids, 2)
    pairs_a <- c(pairs_a, grid[1, ])
    pairs_b <- c(pairs_b, grid[2, ])
    evid <- c(evid, rep(names_vec[i], ncol(grid)))
  }
  pairs <- data.frame(
    uniprot_a = pairs_a,
    uniprot_b = pairs_b,
    evidence = evid,
    stringsAsFactors = FALSE
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
    source_url = paste(
      "https://omnipathdb.org/ (CORUM complexes) /",
      "https://mips.helmholtz-muenchen.de/corum/"
    ),
    license = "CC BY 4.0 (CORUM); see OmniPath source licenses if via OmniPath",
    citation = "CORUM; Tsitsiridis et al. / CORUM release papers"
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
  license <- match.arg(license, choices = c("commercial", "academic", "unknown"))
  if (is.null(interactions_file)) {
    interactions_file <- file.path(
      .default_raw_cache(), "omnipath",
      paste0("interactions_omnipath_", license, ".tsv")
    )
    fallback <- file.path(
      .default_raw_cache(), "omnipath", "interactions_omnipath.tsv"
    )
    if (!file.exists(interactions_file) && file.exists(fallback)) {
      interactions_file <- fallback
      if (license != "unknown") {
        cli::cli_inform(c(
          "i" = "Using cached OmniPath file without re-download: {.path {interactions_file}}",
          "i" = "Requested license filter: {.val {license}}"
        ))
      }
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
  src <- intersect(c("source", "uniprot_a"), colnames(dt))[1]
  tgt <- intersect(c("target", "uniprot_b"), colnames(dt))[1]
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
    source_url = "https://omnipathdb.org/",
    license = paste0("OmniPath per-source; filter=", license),
    citation = "Türei et al.; OmniPath + contributing resources"
  ))
}

#' Max of two directional score columns, ignoring NA on one side
#' @noRd
.max_directional_score <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ifelse(is.na(x), y, ifelse(is.na(y), x, pmax(x, y)))
}

#' Build AlphaFold DB high-confidence complex database
#'
#' When no local CSVs are present, downloads the NVIDIA/AFDB heterodimer
#' metadata table from EBI FTP (~2 GB) into \code{raw_dir}.
#'
#' @param heterodimer_file Optional FTP-derived heterodimer metadata CSV.
#'   If absent, uses panel API cache CSVs under \code{raw_dir}, then downloads
#'   the full heterodimer metadata dump if needed.
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
  files <- unique(c(
    heterodimer_file,
    file.path(raw_dir, "afdb_api_complexes_panel.csv"),
    file.path(raw_dir, "afdb_heterodimers_panel_pairs.csv"),
    file.path(raw_dir, "afdb_homodimers_human_panel.csv"),
    file.path(raw_dir, "heterodimer_metadata.csv")
  ))
  files <- files[file.exists(files)]
  if (length(files) < 1) {
    hetero_dest <- file.path(raw_dir, "heterodimer_metadata.csv")
    cli::cli_inform(c(
      "i" = paste(
        "No local AFDB CSVs found; downloading NVIDIA heterodimer metadata",
        "(large file, ~2 GB)."
      )
    ))
    .download_if_missing(
      paste0(
        "https://ftp.ebi.ac.uk/pub/databases/alphafold/",
        "collaborations/nvda/heterodimer_metadata.csv"
      ),
      hetero_dest
    )
    files <- hetero_dest
  }

  pieces <- lapply(files, function(f) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("tax_id_1", "tax_id_2") %in% names(df))) {
      df <- df[
        as.character(df$tax_id_1) == "9606" &
          as.character(df$tax_id_2) == "9606",
        ,
        drop = FALSE
      ]
    }
    a <- if ("uniprot_ac_1" %in% names(df)) {
      df$uniprot_ac_1
    } else if ("uniprot_a" %in% names(df)) {
      df$uniprot_a
    } else if ("a" %in% names(df)) {
      df$a
    } else {
      NULL
    }
    b <- if ("uniprot_ac_2" %in% names(df)) {
      df$uniprot_ac_2
    } else if ("uniprot_b" %in% names(df)) {
      df$uniprot_b
    } else if ("b" %in% names(df)) {
      df$b
    } else {
      NULL
    }
    if (is.null(a) || is.null(b) || length(a) < 1) {
      return(NULL)
    }
    ipsae <- if ("ipSAE" %in% names(df)) {
      as.numeric(df$ipSAE)
    } else if ("ipsae" %in% names(df)) {
      as.numeric(df$ipsae)
    } else if (all(c("ipSAE_AB", "ipSAE_BA") %in% names(df))) {
      .max_directional_score(df$ipSAE_AB, df$ipSAE_BA)
    } else {
      rep(NA_real_, length(a))
    }
    pdq <- if ("pDockQ2" %in% names(df)) {
      as.numeric(df$pDockQ2)
    } else if ("pdockq2" %in% names(df)) {
      as.numeric(df$pdockq2)
    } else if (all(c("pDockQ2_AB", "pDockQ2_BA") %in% names(df))) {
      .max_directional_score(df$pDockQ2_AB, df$pDockQ2_BA)
    } else {
      rep(NA_real_, length(a))
    }
    data.frame(
      uniprot_a = as.character(a),
      uniprot_b = as.character(b),
      ipSAE = ipsae,
      pDockQ2 = pdq,
      stringsAsFactors = FALSE
    )
  })
  pairs <- bind_rows(pieces)
  if (nrow(pairs) < 1) {
    cli::cli_abort("No usable AFDB pairs found in {.path {raw_dir}}.")
  }
  keep <- (
    (is.na(pairs$ipSAE) & is.na(pairs$pDockQ2)) |
      ((!is.na(pairs$ipSAE) & pairs$ipSAE >= ipsae_min) |
        (!is.na(pairs$pDockQ2) & pairs$pDockQ2 >= pdockq2_min))
  )
  pairs <- pairs[keep, , drop = FALSE]
  pairs$score <- ifelse(!is.na(pairs$ipSAE), pairs$ipSAE, pairs$pDockQ2)
  pairs$evidence <- sprintf("ipSAE=%s;pDockQ2=%s", pairs$ipSAE, pairs$pDockQ2)
  edges <- normalise_interaction_edges(
    pairs,
    score_col = "score",
    evidence_col = "evidence",
    resource = "alphafold",
    resource_version = version
  )

  return(save_interaction_database(
    edges = edges,
    database = "alphafold",
    version = version,
    cache_dir = cache_dir,
    source_url = paste0(
      "https://ftp.ebi.ac.uk/pub/databases/alphafold/",
      "collaborations/nvda/"
    ),
    license = "CC BY 4.0",
    citation = "AlphaFold DB / NVIDIA complexes; Jumper et al. Nature 2021"
  ))
}

#' Build all five interaction databases
#'
#' Maintainer helper that runs each \code{build_*_database()} writer. Missing
#' raw dumps are downloaded into the package raw cache. OmniPath defaults to
#' the commercial license filter.
#'
#' @param cache_dir Output interaction database cache.
#' @return Named list of RDS paths.
#' @export
build_all_interaction_databases <- function(
  cache_dir = interaction_database_cache_dir()
) {
  return(list(
    string = build_string_database(cache_dir = cache_dir),
    biogrid = build_biogrid_database(cache_dir = cache_dir),
    corum = build_corum_database(cache_dir = cache_dir),
    omnipath = build_omnipath_database(
      cache_dir = cache_dir,
      license = "commercial"
    ),
    alphafold = build_alphafold_database(cache_dir = cache_dir)
  ))
}
