#' Build slim interaction priors from public database releases
#'
#' These maintainer/CI helpers download or read local release files and write
#' canonical UniProt edge tables under \code{InteractionPriorCacheDir()}.

.default_raw_cache <- function() {
  if (requireNamespace("here", quietly = TRUE)) {
    return(here::here("results", "db_cache"))
  }
  file.path("results", "db_cache")
}

#' Build STRING physical / full link prior (human)
#'
#' @param raw_dir Directory containing STRING download files (links + aliases).
#' @param version STRING version label (default \code{"12.0"}).
#' @param cache_dir Output prior cache.
#' @param species NCBI taxon (default human 9606).
#' @return Path to saved prior RDS.
#' @export
BuildStringPrior <- function(
  raw_dir = file.path(.default_raw_cache(), "string"),
  version = "12.0",
  cache_dir = InteractionPriorCacheDir(),
  species = 9606,
  include_full = FALSE
) {
  .find_string_file <- function(fname) {
    candidates <- c(
      file.path(raw_dir, fname),
      file.path(.default_raw_cache(), "string", fname),
      file.path(.default_raw_cache(), "alphafold_interactions", fname)
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) hit[[1]] else NA_character_
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
    stop("Missing STRING aliases under ", raw_dir, call. = FALSE)
  }

  aliases <- data.table::fread(
    aliases_gz,
    header = TRUE,
    sep = "\t",
    quote = "",
    showProgress = FALSE
  )
  data.table::setnames(aliases, c("string_protein_id", "alias", "source"))
  # One UniProt accession per STRING id: prefer classical Swiss-Prot ACs
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
      stop("Unexpected STRING links columns in ", path, call. = FALSE)
    }
    m1 <- merge(
      links, up,
      by.x = "protein1", by.y = "string_protein_id", allow.cartesian = TRUE
    )
    data.table::setnames(m1, "uniprot", "uniprot_a")
    m2 <- merge(
      m1, up,
      by.x = "protein2", by.y = "string_protein_id",
      allow.cartesian = TRUE
    )
    data.table::setnames(m2, "uniprot", "uniprot_b")
    edges <- NormaliseInteractionPrior(
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

  # Prefer physical for PNA; full network is large — opt in with include_full=TRUE
  phys <- read_links(physical_gz, "physical")
  full <- if (isTRUE(include_full)) read_links(full_gz, "full") else NULL
  if (is.null(phys) && is.null(full)) {
    stop(
      "No STRING link files found (expected physical and/or full links .txt.gz)",
      call. = FALSE
    )
  }
  edges <- dplyr::bind_rows(phys, full)

  SaveInteractionPrior(
    edges = edges,
    database = "string",
    version = version,
    cache_dir = cache_dir,
    source_url = "https://string-db.org/cgi/download",
    license = "CC BY 4.0",
    citation = "Szklarczyk et al.; STRING database (see string-db.org)"
  )
}

#' Build BioGRID physical interaction prior (human)
#'
#' @param raw_file Path to BioGRID MV-Physical or organism tab3 file.
#' @param version BioGRID release label.
#' @param cache_dir Output prior cache.
#' @return Path to saved prior RDS.
#' @export
BuildBiogridPrior <- function(
  raw_file = file.path(
    .default_raw_cache(), "biogrid",
    "BIOGRID-MV-Physical-5.0.259.tab3.txt"
  ),
  version = "5.0.259",
  cache_dir = InteractionPriorCacheDir()
) {
  if (!file.exists(raw_file)) {
    stop("Missing BioGRID file: ", raw_file, call. = FALSE)
  }
  dt <- data.table::fread(raw_file, sep = "\t", quote = "", showProgress = FALSE)
  # SWISS-PROT accessions
  a_col <- grep("SWISS-PROT.*Interactor A|SWISS-PROT Accessions Interactor A",
    colnames(dt), value = TRUE, ignore.case = TRUE
  )
  b_col <- grep("SWISS-PROT.*Interactor B|SWISS-PROT Accessions Interactor B",
    colnames(dt), value = TRUE, ignore.case = TRUE
  )
  org_a <- grep("Organism ID Interactor A", colnames(dt), value = TRUE)[1]
  org_b <- grep("Organism ID Interactor B", colnames(dt), value = TRUE)[1]
  exp_col <- grep("^Experimental System$", colnames(dt), value = TRUE)[1]
  if (length(a_col) < 1 || length(b_col) < 1) {
    stop("Could not find SWISS-PROT columns in ", raw_file, call. = FALSE)
  }
  a_col <- a_col[1]
  b_col <- b_col[1]

  if (!is.na(org_a) && !is.na(org_b)) {
    dt <- dt[as.character(dt[[org_a]]) == "9606" & as.character(dt[[org_b]]) == "9606"]
  }

  # Take first SWISS-PROT accession per interactor (avoid multi-AC expansion blow-up)
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
  edges <- NormaliseInteractionPrior(
    pairs,
    a_col = "uniprot_a",
    b_col = "uniprot_b",
    evidence_col = "evidence",
    resource = "biogrid",
    resource_version = version
  )

  SaveInteractionPrior(
    edges = edges,
    database = "biogrid",
    version = version,
    cache_dir = cache_dir,
    source_url = "https://downloads.thebiogrid.org/",
    license = "MIT",
    citation = "Stark et al., Nucleic Acids Res. 2006; BioGRID"
  )
}

#' Build CORUM co-membership prior (human)
#'
#' Uses OmniPath-served CORUM complexes when the official CORUM zip is unavailable
#' (previously returned SPA HTML). Prefer an official CORUM file when present.
#'
#' @param corum_file Path to OmniPath CORUM complexes TSV or official coreComplexes.
#' @param version Version label.
#' @param cache_dir Output prior cache.
#' @return Path to saved prior RDS.
#' @export
BuildCorumPrior <- function(
  corum_file = file.path(
    .default_raw_cache(), "corum", "omnipath_complexes_corum.tsv"
  ),
  version = "omnipath_corum",
  cache_dir = InteractionPriorCacheDir()
) {
  if (!file.exists(corum_file)) {
    stop("Missing CORUM complexes file: ", corum_file, call. = FALSE)
  }
  # Guard against SPA HTML stubs
  first <- readLines(corum_file, n = 5, warn = FALSE)
  if (any(grepl("<!DOCTYPE html|<html", first, ignore.case = TRUE))) {
    stop(
      "CORUM file looks like HTML, not data: ", corum_file,
      call. = FALSE
    )
  }

  dt <- data.table::fread(corum_file, sep = "\t", quote = "", showProgress = FALSE)
  # OmniPath complexes: components column with underscore-separated UniProt IDs
  comp_col <- intersect(
    c("components", "subunits_uniprot_id", "subunits"),
    colnames(dt)
  )
  if (length(comp_col) < 1) {
    stop("No components column in ", corum_file, call. = FALSE)
  }
  comp_col <- comp_col[1]
  name_col <- intersect(c("name", "complex_name", "ComplexName"), colnames(dt))
  name_col <- if (length(name_col)) name_col[1] else NULL

  pairs_a <- character()
  pairs_b <- character()
  evid <- character()
  comps <- as.character(dt[[comp_col]])
  names_vec <- if (!is.null(name_col)) as.character(dt[[name_col]]) else rep(NA_character_, length(comps))
  for (i in seq_along(comps)) {
    ids <- unique(trimws(unlist(strsplit(comps[i], "[_;,]"))))
    ids <- ids[ids != "" & !is.na(ids)]
    if (length(ids) < 2) next
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
  edges <- NormaliseInteractionPrior(
    pairs,
    evidence_col = "evidence",
    resource = "corum",
    resource_version = version
  )

  SaveInteractionPrior(
    edges = edges,
    database = "corum",
    version = version,
    cache_dir = cache_dir,
    source_url = "https://omnipathdb.org/ (CORUM complexes) / https://mips.helmholtz-muenchen.de/corum/",
    license = "CC BY 4.0 (CORUM); see OmniPath source licenses if via OmniPath",
    citation = "CORUM; Tsitsiridis et al. / CORUM release papers"
  )
}

#' Build OmniPath interaction prior
#'
#' @param interactions_file Path to OmniPath interactions TSV.
#' @param version Version label.
#' @param cache_dir Output prior cache.
#' @param license OmniPath license filter used when the file was downloaded
#'   (\code{"commercial"} recommended for product builds).
#' @return Path to saved prior RDS.
#' @export
BuildOmnipathPrior <- function(
  interactions_file = NULL,
  version = "webservice",
  cache_dir = InteractionPriorCacheDir(),
  license = c("commercial", "academic", "unknown")
) {
  license <- match.arg(license)
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
        message(
          "Using cached OmniPath file without re-download: ", interactions_file,
          " (requested license filter: ", license, ")"
        )
      }
    }
  }
  if (!file.exists(interactions_file)) {
    if (license == "unknown") {
      stop("Missing OmniPath interactions file: ", interactions_file, call. = FALSE)
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
    interactions_file, sep = "\t", quote = "", showProgress = FALSE
  )
  src <- intersect(c("source", "uniprot_a"), colnames(dt))[1]
  tgt <- intersect(c("target", "uniprot_b"), colnames(dt))[1]
  if (is.na(src) || is.na(tgt)) {
    stop("Could not find source/target UniProt columns", call. = FALSE)
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
  edges <- NormaliseInteractionPrior(
    pairs,
    evidence_col = "evidence",
    resource = "omnipath",
    resource_version = paste0(version, "_", license)
  )

  SaveInteractionPrior(
    edges = edges,
    database = "omnipath",
    version = paste0(version, "_", license),
    cache_dir = cache_dir,
    source_url = "https://omnipathdb.org/",
    license = paste0("OmniPath per-source; filter=", license),
    citation = "Türei et al.; OmniPath + contributing resources"
  )
}

#' Build AlphaFold DB high-confidence complex prior
#'
#' @param heterodimer_file Optional FTP-derived heterodimer metadata CSV
#'   (large). If absent, uses panel API cache CSVs under \code{raw_dir}.
#' @param raw_dir Directory with AFDB cache files.
#' @param version Version label.
#' @param cache_dir Output prior cache.
#' @param ipsae_min Minimum ipSAE.
#' @param pdockq2_min Minimum pDockQ2.
#' @return Path to saved prior RDS.
#' @export
BuildAlphafoldPrior <- function(
  heterodimer_file = file.path(
    .default_raw_cache(), "alphafold_interactions",
    "afdb_heterodimers_human_panel_any.csv"
  ),
  raw_dir = file.path(.default_raw_cache(), "alphafold_interactions"),
  version = "afdb_nvda_highconf",
  cache_dir = InteractionPriorCacheDir(),
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
    stop("No AFDB complex CSVs found under ", raw_dir, call. = FALSE)
  }

  pieces <- lapply(files, function(f) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    a <- df[["uniprot_ac_1"]] %||% df[["uniprot_a"]] %||% df[["a"]]
    b <- df[["uniprot_ac_2"]] %||% df[["uniprot_b"]] %||% df[["b"]]
    if (is.null(a) || is.null(b)) {
      return(NULL)
    }
    ipsae <- df[["ipSAE"]] %||% df[["ipsae"]] %||% NA_real_
    pdq <- df[["pDockQ2"]] %||% df[["pdockq2"]] %||% NA_real_
    data.frame(
      uniprot_a = as.character(a),
      uniprot_b = as.character(b),
      ipSAE = as.numeric(ipsae),
      pDockQ2 = as.numeric(pdq),
      stringsAsFactors = FALSE
    )
  })
  pairs <- dplyr::bind_rows(pieces)
  # keep rows passing either threshold when scores present; keep NA-score rows
  # from sparse exports that were already high-conf filtered
  keep <- (
    (is.na(pairs$ipSAE) & is.na(pairs$pDockQ2)) |
      ((!is.na(pairs$ipSAE) & pairs$ipSAE >= ipsae_min) |
        (!is.na(pairs$pDockQ2) & pairs$pDockQ2 >= pdockq2_min))
  )
  pairs <- pairs[keep, , drop = FALSE]
  pairs$score <- dplyr::coalesce(pairs$ipSAE, pairs$pDockQ2)
  pairs$evidence <- sprintf("ipSAE=%s;pDockQ2=%s", pairs$ipSAE, pairs$pDockQ2)
  edges <- NormaliseInteractionPrior(
    pairs,
    score_col = "score",
    evidence_col = "evidence",
    resource = "alphafold",
    resource_version = version
  )

  SaveInteractionPrior(
    edges = edges,
    database = "alphafold",
    version = version,
    cache_dir = cache_dir,
    source_url = "https://alphafold.ebi.ac.uk/",
    license = "CC BY 4.0",
    citation = "AlphaFold DB / NVIDIA complexes; Jumper et al. Nature 2021"
  )
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Build all five interaction priors from local raw caches
#'
#' @param cache_dir Output prior cache.
#' @return Named list of RDS paths.
#' @export
BuildAllInteractionPriors <- function(cache_dir = InteractionPriorCacheDir()) {
  list(
    string = BuildStringPrior(cache_dir = cache_dir),
    biogrid = BuildBiogridPrior(cache_dir = cache_dir),
    corum = BuildCorumPrior(cache_dir = cache_dir),
    omnipath = BuildOmnipathPrior(cache_dir = cache_dir, license = "commercial"),
    alphafold = BuildAlphafoldPrior(cache_dir = cache_dir)
  )
}
