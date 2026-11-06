#' Read enrichment diagnostics and provenance
#'
#' Metadata is recorded on newly created Annot, richResult, and GSEAResult
#' objects. Legacy objects return an explicit unavailable record.
#' @param x Annot, richResult, or GSEAResult object
#' @details New analyses record schema 2: input occurrences and any matching
#' conversion record, the actual pre-reporting tested family, and SHA-256
#' fingerprints of normalized upstream memberships and resolved term labels.
#' Unknown historical producer or tested evidence is never inferred from results.
#' @return A named list containing diagnostics and provenance
#' @export
enrichmentInfo <- function(x) {
  supported <- methods::is(x, "Annot") || methods::is(x, "richResult") ||
    methods::is(x, "GSEAResult")
  if (!supported) {
    stop("enrichmentInfo() supports Annot, richResult, or GSEAResult objects.",
         call. = FALSE)
  }
  metadata <- attr(x, "richR_metadata", exact = TRUE)
  if (is.null(metadata)) {
    return(list(schema_version = 1L, available = FALSE,
                reason = "metadata not recorded", diagnostics = NULL,
                provenance = NULL))
  }
  metadata
}

.richr_version <- function() {
  if ("richR" %in% loadedNamespaces()) {
    return(as.character(utils::packageVersion("richR")))
  }
  NA_character_
}

.package_version_or_na <- function(package) {
  if (length(package) == 1L && !is.na(package) && nzchar(package) &&
      requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}

.annotation_provenance <- function(object = NULL, constructor = NULL,
                                   species = NULL, keytype = NULL,
                                   anntype = NULL, source_type = NULL,
                                   source_name = NULL, source_version = NA_character_,
                                   provider_package = NULL,
                                   provider_version = NULL,
                                   constructor_normalization = NULL) {
  if (methods::is(object, "Annot")) {
    previous <- attr(object, "richR_metadata", exact = TRUE)
    if (is.list(previous) && is.list(previous$provenance$annotation)) {
      annotation <- previous$provenance$annotation
      annotation$species <- object@species
      annotation$keytype <- object@keytype
      annotation$anntype <- object@anntype
      return(annotation)
    }
    species <- object@species
    keytype <- object@keytype
    anntype <- object@anntype
    source_type <- source_type %||% "unknown"
    source_name <- source_name %||% NA_character_
  }
  list(
    constructor = constructor,
    species = species,
    keytype = keytype,
    anntype = anntype,
    source_type = source_type %||% "custom",
    source_name = source_name %||% NA_character_,
    source_version = source_version,
    provider_package = provider_package,
    provider_version = provider_version,
    constructor_normalization = constructor_normalization
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.annotation_metadata <- function(normalization, provenance, options = list()) {
  list(schema_version = 1L, available = TRUE, analysis = "annotation",
       annotation = normalization$counts,
       provenance = list(richR_version = .richr_version(),
                         R_version = as.character(getRversion()),
                         fgsea_version = NULL,
                         annotation = provenance, rng = NULL),
       options = options)
}

.attach_annotation_metadata <- function(object, normalization, provenance,
                                        options = list()) {
  attr(object, "richR_metadata") <- .annotation_metadata(
    normalization, provenance, options
  )
  object
}

.analysis_metadata <- function(analysis, input, annotation, background,
                               testing, options, annotation_provenance,
                               fgsea_version = NULL, rng = NULL) {
  list(schema_version = 2L, available = TRUE, analysis = analysis,
       input = input, annotation = annotation, background = background,
       testing = testing, options = options,
       provenance = list(richR_version = .richr_version(),
                         R_version = as.character(getRversion()),
                         fgsea_version = fgsea_version,
                         annotation = annotation_provenance,
                         rng = rng))
}


.geneset_fingerprint <- function(annot, term_names = NULL) {
  # Fingerprint actual upstream column semantics; no identity migration.
  frame <- data.frame(GeneID = as.character(annot[[1L]]),
                      id = as.character(annot[[2L]]),
                      label = as.character(annot[[2L]]), stringsAsFactors = FALSE)
  if (!is.null(term_names)) {
    labels <- as.character(unname(term_names[frame$id]))
    missing <- is.na(labels) | !nzchar(trimws(labels))
    labels[missing] <- frame$id[missing]
    frame$label <- labels
  }
  ids <- as.character(frame$id)
  genes <- as.character(frame$GeneID)
  labels <- as.character(frame$label)
  ## Blank = only spaces, TAB, CR or LF (what trimws() strips).
  blank <- function(x) grepl("^[ \t\r\n]*$", x, perl = TRUE)
  valid <- !is.na(ids) & !is.na(genes) & !blank(ids) & !blank(genes)
  ids <- enc2utf8(ids[valid])
  genes <- enc2utf8(genes[valid])
  labels <- labels[valid]
  ## Label of each identifier: its first non-blank label in row order (REQ-003).
  labelled <- which(!is.na(labels) & !blank(labels))
  first <- labelled[!duplicated(ids[labelled])]
  ## Unique pairs in C-locale byte order: sort, then drop repeats of the previous pair.
  o <- order(ids, genes, method = "radix")
  id <- ids[o]
  gene <- genes[o]
  n <- length(id)
  keep <- if (n > 1L) c(TRUE, id[-1L] != id[-n] | gene[-1L] != gene[-n]) else rep(TRUE, n)
  id <- id[keep]
  gene <- gene[keep]
  label <- labels[first][match(id, ids[first])]
  label[is.na(label)] <- id[is.na(label)]
  escape <- function(x) {
    if (!any(grepl("[\\\\\t\n\r]", x, perl = TRUE))) return(x)
    x <- gsub("\\", "\\\\", x, fixed = TRUE)
    x <- gsub("\t", "\\t", x, fixed = TRUE)
    x <- gsub("\n", "\\n", x, fixed = TRUE)
    gsub("\r", "\\r", x, fixed = TRUE)
  }
  id_text <- escape(id)
  gene_text <- escape(gene)
  label <- escape(enc2utf8(label))
  ## Hash the canonical text's bytes, written field by field into memory (much
  ## faster than paste(collapse =) on large annotations; same bytes).
  hash <- function(header, fields) {
    con <- rawConnection(raw(0), "wb")
    on.exit(close(con))
    writeLines(header, con, sep = "\n", useBytes = TRUE)
    if (length(id)) {
      writeLines(as.vector(do.call(rbind, fields)), con, sep = "", useBytes = TRUE)
    }
    digest(rawConnectionValue(con), algo = "sha256", serialize = FALSE)
  }
  m <- length(id)
  list(
    algorithm = "sha256", canonicalization = "richR-geneset-v1",
    membership = hash("richR-geneset-v1", list(id_text, "\t", gene_text, "\n")),
    full = hash("richR-geneset-v1/full",
                list(id_text, "\t", gene_text, "\t", label, "\n")),
    n_pairs = m,
    n_sets = if (m) 1L + sum(id[-1L] != id[-m]) else 0L,
    n_genes = length(unique(gene))
  )
}
