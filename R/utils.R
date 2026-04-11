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
  valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!padj.method %in% valid_methods) {
    stop("padj.method must be one of: ", paste(valid_methods, collapse = ", "))
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
  # Lower bounds: original logic uses OR when keepRich=TRUE
  if (isTRUE(keepRich)) {
    resultFis <- subset(resultFis, Significant >= minSize | RichFactor == 1 | Annotated >= minGSSize)
  } else {
    resultFis <- subset(resultFis, Significant >= minSize)
    resultFis <- subset(resultFis, Annotated >= minGSSize)
  }
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
#' @param annot 2-column data.frame: col1 = term/pathway, col2 = gene
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
#' @return richResult S4 object
#' @keywords internal
.run_ora <- function(x, annot, term_names = NULL, extra_cols = NULL,
                     pvalue = 0.05, padj = NULL, padj.method = "BH",
                     minSize = 2, maxSize = 500, minGSSize = 10, maxGSSize = 500,
                     keepRich = TRUE, organism = NULL, ontology = "",
                     keytype = "", filename = NULL, sep = ",") {
  .validate_enrichment_input(x, annot = annot, pvalue = pvalue, padj = padj,
                             padj.method = padj.method, minSize = minSize, maxSize = maxSize)
  ao2gene <- sf(annot)
  ao2gene_num <- name_table(ao2gene)
  gene2ao <- sf(annot[, c(2, 1)])
  if (is.data.frame(x)) {
    input <- rownames(x)
  } else {
    input <- as.vector(x)
  }
  fgene2ao <- gene2ao[input]
  fao2gene <- reverseList(fgene2ao)
  k <- name_table(fao2gene)
  n <- sum(!is.na(names(fgene2ao)))
  M <- ao2gene_num[names(k)]
  N <- length(unique(unlist(ao2gene)))
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
  .make_richResult(resultFis, detail, pvalue, padj, padj.method,
                   input, organism, ontology, keytype, sep)
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
