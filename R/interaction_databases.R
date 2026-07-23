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
  edges <- edges[
    edges$uniprot_a %in% panel_uniprot & edges$uniprot_b %in% panel_uniprot,
    ,
    drop = FALSE
  ]
  if (nrow(edges) == 0) {
    return(.empty_panel_interactions())
  }

  m_a <- merge(
    edges,
    map,
    by.x = "uniprot_a",
    by.y = "uniprot_id"
  )
  names(m_a)[names(m_a) == "marker"] <- "marker_1"
  m_ab <- merge(
    m_a,
    map,
    by.x = "uniprot_b",
    by.y = "uniprot_id"
  )
  names(m_ab)[names(m_ab) == "marker"] <- "marker_2"

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
#' Maintainer helper: reads local STRING release files and writes a slim
#' UniProt edge cache via \code{\link{save_interaction_database}}.
#'
#' @param raw_dir Directory containing STRING download files (links + aliases).
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

  .find_string_file <- function(fname) {
    path <- file.path(raw_dir, fname)
    if (file.exists(path)) path else NA_character_
  }

  aliases_gz <- .find_string_file(
    sprintf("%s.protein.aliases.v%s.txt.gz", species, version)
  )
  physical_gz <- .find_string_file(
    sprintf("%s.protein.physical.links.v%s.txt.gz", species, version)
  )
  full_gz <- .find_string_file(
    sprintf("%s.protein.links.v%s.txt.gz", species, version)
  )

  if (is.na(aliases_gz)) {
    cli::cli_abort("Missing STRING aliases under {.path {raw_dir}}.")
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
    if (is.na(path) || !file.exists(path)) {
      return(NULL)
    }
    links <- data.table::fread(path, header = TRUE, sep = " ", showProgress = FALSE)
    data.table::setnames(links, tolower(gsub("#", "", names(links))))
    if (!all(c("protein1", "protein2", "combined_score") %in% names(links))) {
      cli::cli_abort("Unexpected STRING links columns in {.path {path}}.")
    }
    m1 <- merge(
      links, up,
      by.x = "protein1", by.y = "string_protein_id"
    )
    data.table::setnames(m1, "uniprot", "uniprot_a")
    m2 <- merge(
      m1, up,
      by.x = "protein2", by.y = "string_protein_id"
    )
    data.table::setnames(m2, "uniprot", "uniprot_b")
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
#' @param raw_file Path to BioGRID MV-Physical or organism tab3 file.
#' @param version BioGRID release label.
#' @param cache_dir Output interaction database cache.
#' @return Path to saved RDS.
#' @export
build_biogrid_database <- function(
  raw_file = file.path(
    .default_raw_cache(), "biogrid",
    "BIOGRID-MV-Physical-5.0.259.tab3.txt"
  ),
  version = "5.0.259",
  cache_dir = interaction_database_cache_dir()
) {
  rlang::check_installed("data.table")
  assert_file_exists(raw_file)
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
#' Uses OmniPath-served CORUM complexes when the official CORUM zip is
#' unavailable. Prefer an official CORUM file when present.
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
  assert_file_exists(corum_file)
  first <- readLines(corum_file, n = 5, warn = FALSE)
  if (any(grepl("<!DOCTYPE html|<html", first, ignore.case = TRUE))) {
    cli::cli_abort("CORUM file looks like HTML, not data: {.path {corum_file}}.")
  }

  dt <- data.table::fread(corum_file, sep = "\t", quote = "", showProgress = FALSE)
  comp_col <- intersect(
    c("components", "subunits_uniprot_id", "subunits"),
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
    dir.create(dirname(interactions_file), recursive = TRUE, showWarnings = FALSE)
    utils::download.file(url, interactions_file, mode = "wb", quiet = TRUE)
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

#' Build AlphaFold DB high-confidence complex database
#'
#' @param heterodimer_file Optional FTP-derived heterodimer metadata CSV.
#'   If absent, uses panel API cache CSVs under \code{raw_dir}.
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
    file.path(raw_dir, "afdb_homodimers_human_panel.csv")
  ))
  files <- files[file.exists(files)]
  if (length(files) < 1) {
    cli::cli_abort("No AFDB complex CSVs found under {.path {raw_dir}}.")
  }

  pieces <- lapply(files, function(f) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    a <- if ("uniprot_ac_1" %in% names(df)) {
      df$uniprot_ac_1
    } else if ("uniprot_a" %in% names(df)) {
      df$uniprot_a
    } else {
      df$a
    }
    b <- if ("uniprot_ac_2" %in% names(df)) {
      df$uniprot_ac_2
    } else if ("uniprot_b" %in% names(df)) {
      df$uniprot_b
    } else {
      df$b
    }
    if (is.null(a) || is.null(b)) {
      return(NULL)
    }
    ipsae <- if ("ipSAE" %in% names(df)) {
      df$ipSAE
    } else if ("ipsae" %in% names(df)) {
      df$ipsae
    } else {
      NA_real_
    }
    pdq <- if ("pDockQ2" %in% names(df)) {
      df$pDockQ2
    } else if ("pdockq2" %in% names(df)) {
      df$pdockq2
    } else {
      NA_real_
    }
    data.frame(
      uniprot_a = as.character(a),
      uniprot_b = as.character(b),
      ipSAE = as.numeric(ipsae),
      pDockQ2 = as.numeric(pdq),
      stringsAsFactors = FALSE
    )
  })
  pairs <- bind_rows(pieces)
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
    source_url = "https://alphafold.ebi.ac.uk/",
    license = "CC BY 4.0",
    citation = "AlphaFold DB / NVIDIA complexes; Jumper et al. Nature 2021"
  ))
}

#' Build all five interaction databases from local raw caches
#'
#' Maintainer helper that runs each \code{build_*_database()} writer. OmniPath
#' defaults to the commercial license filter.
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
