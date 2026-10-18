## ----------------------------------------------------------------
##  ORA helper functions (shared across enrichment functions)
## ----------------------------------------------------------------

#' Validate common enrichment inputs
#' @param x gene vector or data.frame
#' @param annot annotation data.frame (2+ columns: GeneID, Term)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff (or NULL)
#' @param padj.method p-value adjustment method
#' @param minSize minimum significant genes
#' @param maxSize maximum significant genes
#' @keywords internal
.validate_enrichment_input <- function(x, annot = NULL, pvalue = 0.05, padj = NULL,
                                        padj.method = "BH", minSize = 1, maxSize = 500) {
  # Validate gene input
  if (is.data.frame(x)) {
    genes <- rownames(x)
  } else {
    genes <- as.vector(x)
  }
  if (length(genes) == 0) {
    stop("Input gene vector is empty")
  }
  if (all(is.na(genes))) {
    stop("Input gene vector contains only NA values")
  }
  if (any(duplicated(genes[!is.na(genes)]))) {
    warning("Duplicated genes detected in input; duplicates will be removed")
  }
  # Validate p-value parameters
  if (!is.null(pvalue) && (pvalue < 0 || pvalue > 1)) {
    stop("pvalue must be between 0 and 1")
  }
  if (!is.null(padj) && (padj < 0 || padj > 1)) {
    stop("padj must be between 0 and 1")
  }
  # Validate adjustment method
  if (length(padj.method) != 1L || is.na(padj.method) ||
      !padj.method %in% stats::p.adjust.methods) {
    stop("padj.method must be one of: ", paste(stats::p.adjust.methods, collapse = ", "))
  }
  # Validate size
  if (minSize > maxSize) {
    stop("minSize must be <= maxSize")
  }
  # Validate annotation
  if (!is.null(annot)) {
    if (is.data.frame(annot) && ncol(annot) < 2) {
      stop("Annotation data.frame must have at least 2 columns")
    }
    if (is.data.frame(annot) && nrow(annot) == 0) {
      stop("Annotation data.frame is empty")
    }
  }
  invisible(TRUE)
}

#' Normalize annotation membership rows for enrichment
#' @param annot annotation data.frame with gene and term in its first two columns
#' @param func_name name used in validation messages
#' @return normalized data.frame and row counts
#' @keywords internal
.normalize_annotation <- function(annot, func_name = "enrichment") {
  if (!is.data.frame(annot)) {
    stop(func_name, ": annotation must be a data.frame.", call. = FALSE)
  }
  if (ncol(annot) < 2L) {
    stop(func_name, ": annotation must have at least 2 columns (gene and term).",
         call. = FALSE)
  }
  if (nrow(annot) == 0L) {
    stop(func_name, ": annotation has no membership rows.", call. = FALSE)
  }
  valid_col <- vapply(annot[1:2], function(col) {
    is.atomic(col) && is.null(dim(col))
  }, logical(1))
  if (!all(valid_col)) {
    stop(func_name, ": annotation gene and term columns must be atomic vectors.",
         call. = FALSE)
  }

  genes <- as.character(annot[[1]])
  terms <- as.character(annot[[2]])
  valid <- !is.na(genes) & !is.na(terms) &
    nzchar(trimws(genes)) & nzchar(trimws(terms))
  original_rows <- nrow(annot)
  valid_rows <- sum(valid)
  if (valid_rows == 0L) {
    stop(func_name, ": annotation has no valid annotation memberships.", call. = FALSE)
  }
  invalid_rows <- original_rows - valid_rows
  if (invalid_rows > 0L) {
    warning(func_name, ": excluded ", invalid_rows,
            " invalid annotation row(s) with missing or blank gene/term identifiers.",
            call. = FALSE)
  }

  normalized <- annot[valid, , drop = FALSE]
  normalized[[1]] <- genes[valid]
  normalized[[2]] <- terms[valid]
  duplicated_pair <- duplicated(data.frame(
    gene = normalized[[1]], term = normalized[[2]], stringsAsFactors = FALSE
  ))
  duplicate_rows <- sum(duplicated_pair)
  normalized <- normalized[!duplicated_pair, , drop = FALSE]
  rownames(normalized) <- NULL
  list(
    data = normalized,
    counts = list(
      original_rows = as.integer(original_rows),
      valid_rows = as.integer(valid_rows),
      invalid_rows = as.integer(invalid_rows),
      duplicate_rows = as.integer(duplicate_rows),
      unique_rows = as.integer(nrow(normalized))
    )
  )
}

#' Validate exact identities on a named list before processing its members
#' @param x named list
#' @param arg argument name for messages
#' @return invisible TRUE when names are valid
#' @keywords internal
.validate_named_list <- function(x, arg) {
  if (!is.list(x) || length(x) == 0L || is.null(names(x))) {
    stop(arg, " must be a non-empty named list.", call. = FALSE)
  }
  ids <- names(x)
  valid_name <- !is.na(ids)
  valid_name[valid_name] <- nzchar(ids[valid_name]) & nzchar(trimws(ids[valid_name]))
  invalid <- !valid_name
  if (any(invalid)) {
    stop(arg, " has missing or blank names at position ",
         paste(which(invalid), collapse = ", "), call. = FALSE)
  }
  duplicate <- duplicated(ids) | duplicated(ids, fromLast = TRUE)
  if (any(duplicate)) {
    stop(arg, " has duplicate name(s) at position(s) ",
         paste(which(duplicate), collapse = ", "), ": ",
         paste(unique(ids[duplicate]), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Compute ORA statistics (shared across enrichment functions)
#' @param k named numeric vector of significant gene counts per term
#' @param M named numeric vector of annotated gene counts per term
#' @param N total number of annotated genes (background)
#' @param n number of input genes found in annotation
#' @param padj.method p-value adjustment method (default: "BH")
#' @return list with rhs (p-values), lhs (adjusted p-values), Annotated, Significant,
#'   RichFactor, FoldEnrichment, zscore
#' @keywords internal
.compute_ora_stats <- function(k, M, N, n, padj.method = "BH") {
  # k = significant genes per term, M = annotated genes per term
  # N = total annotated genes, n = total input genes in annotation
  rhs <- hyper_bench_vector(k, M, N, n)
  lhs <- stats::p.adjust(rhs, method = padj.method)
  RichFactor <- k / M
  FoldEnrichment <- RichFactor * N / n
  # z-score (hypergeometric distribution mean and variance)
  mu <- M * n / N
  sigma <- mu * (N - n) * (N - M) / N / (N - 1)
  zscore <- (k - mu) / sqrt(sigma)
  zscore[!is.finite(zscore)] <- NA_real_
  list(rhs = rhs, lhs = lhs, Annotated = M, Significant = k,
       RichFactor = RichFactor, FoldEnrichment = FoldEnrichment, zscore = zscore)
}

#' Filter ORA result data.frame (shared across enrichment functions)
#' @param resultFis data.frame with enrichment results
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff (or NULL)
#' @param minSize minimum significant genes
#' @param maxSize maximum significant genes
#' @param minGSSize minimum annotated genes for testing
#' @param maxGSSize maximum annotated genes for testing
#' @param keepRich keep terms with RichFactor == 1
#' @return filtered data.frame
#' @keywords internal
.filter_ora_result <- function(resultFis, pvalue, padj, minSize, maxSize,
                                minGSSize, maxGSSize, keepRich) {
  if (nrow(resultFis) == 0) return(resultFis)
  resultFis <- resultFis[order(resultFis$Pvalue), ]
  # Upper bounds always apply
  resultFis <- subset(resultFis, Significant <= maxSize)
  resultFis <- subset(resultFis, Annotated <= maxGSSize)
  lower_bounds <- (resultFis$Significant >= minSize & resultFis$Annotated >= minGSSize) |
    (isTRUE(keepRich) & resultFis$RichFactor == 1)
  resultFis <- resultFis[lower_bounds, , drop = FALSE]
  # Filter by p-value
  if (!is.null(padj)) {
    resultFis <- resultFis[resultFis$Padj < padj, ]
  } else {
    resultFis <- resultFis[resultFis$Pvalue < pvalue, ]
  }
  rownames(resultFis) <- resultFis$Annot
  resultFis
}

#' Build detail data.frame from enrichment results (shared across enrichment functions)
#' @param resultFis enrichment result data.frame
#' @param x original input (gene vector or data.frame)
#' @param sep separator for gene IDs
#' @return data.frame with gene-term detail mapping
#' @keywords internal
.build_detail <- function(resultFis, x, sep = ",") {
  if (nrow(resultFis) == 0) {
    return(data.frame(TERM = character(), Annot = character(),
                      GeneID = character(), Pvalue = numeric(),
                      Padj = numeric(), stringsAsFactors = FALSE))
  }
  if (is.data.frame(x)) {
    detail <- getdetail(resultFis, x, sep = sep)
  } else {
    gene_lists <- strsplit(as.character(resultFis$GeneID), split = sep)
    names(gene_lists) <- resultFis$Annot
    detail <- data.frame(
      TERM = rep(names(gene_lists), times = vapply(gene_lists, length, integer(1))),
      Annot = rep(resultFis$Term[match(names(gene_lists), resultFis$Annot)],
                  times = vapply(gene_lists, length, integer(1))),
      GeneID = unlist(gene_lists),
      Pvalue = rep(resultFis$Pvalue, times = vapply(gene_lists, length, integer(1))),
      Padj = rep(resultFis$Padj, times = vapply(gene_lists, length, integer(1))),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    detail$GeneID <- as.character(detail$GeneID)
  }
  detail
}

#' Create a richResult S4 object (shared across enrichment functions)
#' @param resultFis enrichment result data.frame
#' @param detail detail data.frame
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param padj.method p-value adjustment method
#' @param input character vector of input genes
#' @param organism organism name
#' @param ontology ontology type
#' @param keytype key type
#' @param sep separator
#' @return richResult S4 object
#' @keywords internal
.make_richResult <- function(resultFis, detail, pvalue, padj, padj.method,
                              input, organism, ontology, keytype, sep = ",") {
  if (is.null(organism)) organism <- character()
  if (is.null(ontology)) ontology <- character()
  if (is.null(keytype)) keytype <- character()
  if (is.null(padj)) padj <- numeric()
  new("richResult",
      result = resultFis,
      detail = detail,
      pvalueCutoff = pvalue,
      pAdjustMethod = padj.method,
      padjCutoff = padj,
      genenumber = length(input),
      organism = organism,
      ontology = ontology,
      gene = input,
      keytype = keytype,
      sep = sep
  )
}

#' Run the shared ORA enrichment pipeline
#'
#' Core ORA logic shared by enrich_internal, richGO_internal, richKEGG_internal.
#' Handles: input extraction, annotation mapping via sf/reverseList, hypergeometric
#' test, result assembly, filtering, detail building, and richResult construction.
#'
#' @param x gene vector or data.frame (rownames used if data.frame)
#' @param annot data.frame with gene in column 1 and term/pathway in column 2
#' @param term_names named vector mapping term IDs to readable names (or NULL)
#' @param extra_cols data.frame to cbind to result before filtering (e.g. KEGG path levels)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff (or NULL)
#' @param padj.method p-value adjustment method
#' @param minSize minimum significant genes per term
#' @param maxSize maximum significant genes per term
#' @param minGSSize minimum annotated genes for testing
#' @param maxGSSize maximum annotated genes for testing
#' @param keepRich keep terms with RichFactor == 1
#' @param organism organism name
#' @param ontology ontology type
#' @param keytype key type
#' @param filename optional output filename
#' @param sep separator for gene IDs
#' @param universe optional character vector defining the ORA background
#' @param annotation_provenance optional annotation provenance captured from an Annot object
#' @param metadata_ontology ontology label retained for diagnostics
#' @return richResult S4 object
#' @keywords internal
.run_ora <- function(x, annot, term_names = NULL, extra_cols = NULL,
                     pvalue = 0.05, padj = NULL, padj.method = "BH",
                     minSize = 2, maxSize = 500, minGSSize = 10, maxGSSize = 500,
                     keepRich = TRUE, organism = NULL, ontology = "",
                     keytype = "", filename = NULL, sep = ",", universe = NULL,
                     annotation_provenance = NULL, metadata_ontology = ontology) {
  normalization <- if (!is.data.frame(annot) && is.list(annot) &&
                       all(c("data", "counts") %in% names(annot))) {
    annot
  } else {
    .normalize_annotation(annot, func_name = "ORA")
  }
  annot <- normalization$data
  .validate_enrichment_input(x, annot = annot, pvalue = pvalue, padj = padj,
                             padj.method = padj.method, minSize = minSize, maxSize = maxSize)
  raw_input <- if (is.data.frame(x)) rownames(x) else as.vector(x)
  raw_input <- as.character(raw_input)
  valid_input <- raw_input[!is.na(raw_input) & nzchar(trimws(raw_input))]
  input <- unique(valid_input)
  annotation_genes <- unique(as.character(annot[, 1]))
  annotation_counts <- normalization$counts
  if (is.null(annotation_counts$eligible_rows)) {
    annotation_counts$eligible_rows <- as.integer(nrow(annot))
  }
  if (is.null(universe)) {
    effective_universe <- annotation_genes
    background_counts <- list(
      mode = "annotation", original_count = NA_integer_, valid_count = NA_integer_,
      unique_count = NA_integer_, invalid_count = NA_integer_,
      duplicate_count = NA_integer_, unmapped_count = NA_integer_,
      unmapped_ids = character(), effective_count = as.integer(length(effective_universe)),
      effective_ids = effective_universe
    )
  } else {
    if (!is.atomic(universe) || !is.character(universe)) {
      stop("ORA: universe must be an atomic character vector.", call. = FALSE)
    }
    if (!is.null(dim(universe))) {
      stop("ORA: universe must not have dimensions.", call. = FALSE)
    }
    if (length(universe) == 0) {
      stop("ORA: universe is explicitly empty.", call. = FALSE)
    }
    clean_universe <- unique(universe[!is.na(universe) & nzchar(trimws(universe))])
    if (length(clean_universe) == 0) {
      stop("ORA: universe contains no valid genes after removing NA and blank values.", call. = FALSE)
    }
    effective_universe <- clean_universe[clean_universe %in% annotation_genes]
    if (length(effective_universe) == 0) {
      stop("ORA: universe has no annotated genes after intersecting with the eligible annotation.", call. = FALSE)
    }
    valid_universe <- universe[!is.na(universe) & nzchar(trimws(universe))]
    background_counts <- list(
      mode = "explicit", original_count = as.integer(length(universe)),
      valid_count = as.integer(length(valid_universe)),
      unique_count = as.integer(length(clean_universe)),
      invalid_count = as.integer(length(universe) - length(valid_universe)),
      duplicate_count = as.integer(length(valid_universe) - length(clean_universe)),
      unmapped_count = as.integer(sum(!clean_universe %in% annotation_genes)),
      unmapped_ids = clean_universe[!clean_universe %in% annotation_genes],
      effective_count = as.integer(length(effective_universe)),
      effective_ids = effective_universe
    )
  }
  effective_input <- input[input %in% effective_universe]
  if (length(effective_input) == 0) {
    stop("ORA: no query genes remain in the effective universe.", call. = FALSE)
  }
  annot <- annot[as.character(annot[, 1]) %in% effective_universe, , drop = FALSE]
  annotation_counts$effective_rows <- as.integer(nrow(annot))
  ao2gene <- sf(annot)
  ao2gene_num <- name_table(ao2gene)
  gene2ao <- sf(annot[, c(2, 1)])
  fgene2ao <- gene2ao[effective_input]
  fao2gene <- reverseList(fgene2ao)
  k <- name_table(fao2gene)
  n <- length(effective_input)
  M <- ao2gene_num[names(k)]
  N <- length(effective_universe)
  stats <- .compute_ora_stats(k, M, N, n, padj.method)
  rhs_gene <- unlist(lapply(fao2gene, function(g) paste(unique(g), sep = "", collapse = sep)))
  GeneID <- rhs_gene[names(stats$rhs)]
  # Term names: use supplied mapping or fall back to term IDs
  if (!is.null(term_names)) {
    Term <- term_names[names(stats$rhs)]
  } else {
    Term <- names(stats$rhs)
  }
  resultFis <- data.frame(
    Annot = names(stats$rhs), Term = Term,
    Annotated = stats$Annotated, Significant = stats$Significant,
    RichFactor = stats$RichFactor, FoldEnrichment = stats$FoldEnrichment,
    zscore = stats$zscore, Pvalue = as.vector(stats$rhs), Padj = stats$lhs,
    GeneID = GeneID
  )
  if (!is.null(extra_cols)) {
    resultFis <- cbind(resultFis, extra_cols[resultFis$Annot, , drop = FALSE])
  }
  resultFis <- .filter_ora_result(resultFis, pvalue, padj, minSize, maxSize,
                                   minGSSize, maxGSSize, keepRich)
  if (!is.null(filename)) {
    write.table(resultFis, file = paste(filename, ".txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  detail <- .build_detail(resultFis, x, sep)
  result <- .make_richResult(resultFis, detail, pvalue, padj, padj.method,
                             input, organism, ontology, keytype, sep)
  input_metadata <- list(
    original_count = as.integer(length(raw_input)),
    valid_count = as.integer(length(valid_input)),
    unique_count = as.integer(length(input)),
    duplicate_count = as.integer(length(valid_input) - length(input)),
    invalid_count = as.integer(length(raw_input) - length(valid_input)),
    ids = input,
    mapped_ids = input[input %in% annotation_genes],
    unmapped_ids = input[!input %in% annotation_genes],
    mapped_count = as.integer(sum(input %in% annotation_genes)),
    unmapped_count = as.integer(sum(!input %in% annotation_genes)),
    effective_ids = effective_input,
    effective_count = as.integer(length(effective_input)),
    excluded_by_universe_ids = setdiff(input[input %in% annotation_genes], effective_universe)
  )
  testing <- list(family_policy = "positive_overlap_before_reporting_filters",
                  tested_count = as.integer(length(stats$rhs)),
                  effective_N = as.integer(N), effective_n = as.integer(n))
  options <- list(pvalue = pvalue, padj = padj, minSize = minSize,
                  maxSize = maxSize, minGSSize = minGSSize,
                  maxGSSize = maxGSSize, keepRich = keepRich,
                  padj.method = padj.method, ontology = metadata_ontology,
                  keytype = keytype, universe = if (is.null(universe)) NULL else universe)
  if (is.null(annotation_provenance)) {
    annotation_provenance <- .annotation_provenance(
      species = organism, keytype = keytype, anntype = metadata_ontology,
      source_type = "custom", source_name = NA_character_
    )
  }
  attr(result, "richR_metadata") <- .analysis_metadata(
    "ORA", input_metadata, annotation_counts, background_counts, testing,
    options, annotation_provenance
  )
  result
}

## ----------------------------------------------------------------
##  String utility functions
## ----------------------------------------------------------------

#' Insert newlines after every n words (defaults to 4)
#'
#' @param x Character vector to process.
#' @param n Integer. Number of words per line (default 4).
#' @return Character vector with inserted newlines.
#' @keywords internal
.paste.char <- function(x, n = 4) {
  sapply(x, function(elem) {
    words <- strsplit(elem, "\\s+")[[1]]
    if (length(words) <= n) {
      return(elem)
    }
    groups <- split(words, ceiling(seq_along(words) / n))
    grouped_strings <- sapply(groups, paste, collapse = " ")
    paste(grouped_strings, collapse = "\n")
  }, USE.NAMES = FALSE)
}

#' Convert a named vector to a two-column data.frame
#' @param x named vector
#' @param name character vector of length 2 for column names
#' @return data.frame
#' @keywords internal
vec_to_df <- function(x, name) {
  dd <- data.frame(names(x), x, stringsAsFactors = FALSE)
  colnames(dd) <- name
  dd
}

#' Color scale helper (modified from DOSE)
#' @param c1 start color
#' @param c2 end color
#' @return character vector of 200 colors
#' @keywords internal
.color_scale <- function(c1 = "pink", c2 = "red") {
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(200)
  return(colors)
}

#' Get index for color mapping (modified from DOSE)
#' @param v numeric value
#' @param MIN minimum value
#' @param MAX maximum value
#' @return integer index
#' @keywords internal
.getIdx <- function(v, MIN, MAX) {
  intervals <- seq(MIN, MAX, length.out = 200)
  max(which(intervals <= v))
}

#' Jaccard overlap between two vectors
#' @param x first vector or list
#' @param y second vector or list
#' @return numeric Jaccard index
#' @keywords internal
overlap <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y)) / length(unique(c(x, y)))
}

#' Get GO annotation terms
#' @importFrom AnnotationDbi keys
#' @param ont GO ontology: "BP", "CC", or "MF"
#' @return data.frame with rownames as GO IDs and "annotation" column
#' @keywords internal
.get_go_dat <- function(ont = "BP") {
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("GO.db package required")
  godb <- getExportedValue("GO.db", "GO.db")
  key <- keys(godb)
  suppressMessages(go_dat <- AnnotationDbi::select(godb, keys = key,
                                                    columns = c("TERM", "ONTOLOGY"),
                                                    keytype = "GOID"))
  if (ont == "BP") res <- as.data.frame(subset(go_dat, ONTOLOGY == "BP"))
  if (ont == "CC") res <- as.data.frame(subset(go_dat, ONTOLOGY == "CC"))
  if (ont == "MF") res <- as.data.frame(subset(go_dat, ONTOLOGY == "MF"))
  rownames(res) <- res[, 1]
  res <- res[, 2, drop = FALSE]
  colnames(res) <- "annotation"
  return(res)
}

#' Get KEGG pathway annotation terms
#' @importFrom KEGGREST keggList
#' @param builtin use built-in KEGG data or fetch from API
#' @return data.frame with rownames as pathway IDs and "annotation" column
#' @keywords internal
.get_kg_dat <- function(builtin = TRUE) {
  if (isTRUE(builtin)) {
    kegg_env <- new.env(parent = emptyenv())
    data("kegg", envir = kegg_env)
    return(kegg_env$kegg.db)
  } else {
    pathway <- cbind(keggList('pathway'))
    rownames(pathway) <- sub('.*map', '', rownames(pathway))
    colnames(pathway) <- "annotation"
    pathway <- as.data.frame(pathway)
    pathway$annotation <- as.vector(pathway$annotation)
    if (is.na(pathway['04148', ])) {
      kegg_env <- new.env(parent = emptyenv())
      data("kegg", envir = kegg_env)
      pathway <- kegg_env$kegg.db
    }
    return(pathway)
  }
}

#' Get KEGG module annotation terms
#' @importFrom KEGGREST keggList
#' @return data.frame module annotation
#' @keywords internal
.get_kgm.data <- function() {
  mod_env <- new.env(parent = emptyenv())
  data("module", envir = mod_env)
  return(mod_env$module)
}

#' Build annotation lookup for GO, KEGG, or Module
#' @param ontype one of "GO", "KEGG", "Module"
#' @param builtin use built-in annotation or not (default TRUE)
#' @return data.frame with annotation terms
#' @export
#' @author Kai Guo
getann <- function(ontype = "GO", builtin = TRUE) {
  if (ontype == "GO") {
    res <- rbind(.get_go_dat("BP"), .get_go_dat("MF"), .get_go_dat("CC"))
  }
  if (ontype == "KEGG") {
    res <- .get_kg_dat(builtin = builtin)
  }
  if (ontype == "Module") {
    res <- .get_kgm.data()
  }
  return(res)
}

#' Get all children terms of a GO node
#' @param node GO term ID (default: "GO:0008150" for BP root)
#' @param ontology one of "BP", "CC", "MF"
#' @return character vector of child GO IDs
#' @author Kai Guo
GO_child <- function(node = "GO:0008150", ontology = "BP") {
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("GO.db package required")
  if (ontology == "BP") res <- c(node, getExportedValue("GO.db", "GOBPOFFSPRING")[[node]])
  if (ontology == "CC") res <- c(node, getExportedValue("GO.db", "GOCCOFFSPRING")[[node]])
  if (ontology == "MF") res <- c(node, getExportedValue("GO.db", "GOMFOFFSPRING")[[node]])
  return(res[!is.na(res)])
}

#' Get detail and integrate with the input gene information
#' @importFrom dplyr left_join
#' @param rese richResult or GSEAResult
#' @param resd dataframe with input gene as rownames
#' @param sep character string used to separate the genes when concatenating
#' @return data.frame with gene-term detail mapping
#' @examples
#' \dontrun{
#'   hsako <- buildAnnot(species="human", keytype="SYMBOL", anntype="KEGG")
#'   hsako <- as.data.frame(hsako)
#'   gene <- sample(unique(hsako$GeneID), 1000)
#'   res <- richKEGG(gene, kodata=hsako)
#'   gened <- data.frame(lfc=rnorm(length(gene)))
#'   rownames(gened) <- gene
#'   head(getdetail(res, gened, sep=","))
#' }
#' @export
#' @author Kai Guo
getdetail <- function(rese, resd, sep) {
  if (!is.data.frame(resd)) {
    resd <- data.frame(gene = resd)
  }
  if (!("gene" %in% colnames(resd))) {
    resd$gene <- rownames(resd)
  }
  gene <- strsplit(as.vector(rese$GeneID), split = sep)
  names(gene) <- rese$Annot
  gened <- data.frame(
    "TERM" = rep(names(gene), times = unlist(lapply(gene, length))),
    "Annot" = rep(rese$Term, times = unlist(lapply(gene, length))),
    "GeneID" = unlist(gene), row.names = NULL,
    "Pvalue" = rep(rese$Pvalue, times = unlist(lapply(gene, length))),
    "Padj" = rep(rese$Padj, times = unlist(lapply(gene, length)))
  )
  gened$GeneID <- as.character(gened$GeneID)
  res <- left_join(gened, resd, by = c("GeneID" = "gene"))
  return(res)
}

#' Print MSIGDB category information
#' @export
msigdbinfo <- function() {
  cat("#--------------------------------------------------------------#\n")
  cat("# Molecular Signatures Database                                #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Category | Subcategory # Details ----------------------------#\n")
  cat("# C1               # Positional                                #\n")
  cat("# C2 | CGP         # Chemical and Genetic Perturbations        #\n")
  cat("# C2 | CP          # Canonical Pathways                        #\n")
  cat("# C2 | BIOCARTA    # Canonical BIOCARTA                        #\n")
  cat("# C2 | KEGG        # Canonical KEGG                            #\n")
  cat("# C2 | REACTOME    # Canonical REACTOME                        #\n")
  cat("# C3 | MIR         # Motif miRNA Targets                       #\n")
  cat("# C3 | TFT         # Motif Transcription Factor Targets        #\n")
  cat("# C4 | CGN         # Cancer Gene Neighborhoods                 #\n")
  cat("# C4 | CM          # Cancer Modules                            #\n")
  cat("# C5 | BP          # GO Biological Process                     #\n")
  cat("# C5 | CC          # GO Cellular Component                     #\n")
  cat("# C5 | MF          # GO Molecular Function                     #\n")
  cat("# C6               # Oncogenic Signatures                      #\n")
  cat("# C7               # Immunologic Signatures                    #\n")
  cat("# H                # Hallmark                                   #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Source: http://software.broadinstitute.org/gsea/msigdb       #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Supported species (via msigdbr):                             #\n")
  tbl <- .species_table()
  spe <- tbl$species[!is.na(tbl$msigdb)]
  cat(sort(spe), "\n")
}
