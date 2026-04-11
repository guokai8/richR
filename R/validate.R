#' Validate gene input for enrichment analysis
#'
#' Checks that gene input is non-empty, removes NAs and duplicates,
#' and reports mapping statistics.
#'
#' @param x gene vector or data.frame with gene rownames
#' @param annotation annotation data.frame (two-column: GeneID, Term)
#' @param func_name name of the calling function (for messages)
#' @param verbose logical, print mapping statistics
#' @return cleaned character vector of gene IDs
#' @keywords internal
.validateGeneInput <- function(x, annotation = NULL, func_name = "enrich", verbose = TRUE) {
  if (is.data.frame(x)) {
    input <- rownames(x)
  } else {
    input <- as.vector(x)
  }
  if (length(input) == 0 || all(is.na(input))) {
    stop(func_name, ": input gene list is empty or all NA.", call. = FALSE)
  }
  # Remove NAs
  n_na <- sum(is.na(input))
  input <- input[!is.na(input)]
  # Remove empty strings
  input <- input[nchar(input) > 0]
  # Remove duplicates
  n_before <- length(input)
  input <- unique(input)
  n_dup <- n_before - length(input)
  if (length(input) == 0) {
    stop(func_name, ": no valid genes remaining after removing NAs and empty values.", call. = FALSE)
  }
  if (verbose) {
    if (n_na > 0) message(func_name, ": removed ", n_na, " NA values from input.")
    if (n_dup > 0) message(func_name, ": removed ", n_dup, " duplicate genes.")
  }
  # Check overlap with annotation
  if (!is.null(annotation)) {
    annot_genes <- unique(as.character(annotation[, 1]))
    n_mapped <- sum(input %in% annot_genes)
    pct <- round(100 * n_mapped / length(input), 1)
    if (verbose) {
      message(func_name, ": ", n_mapped, "/", length(input),
              " (", pct, "%) input genes found in annotation.")
    }
    if (n_mapped == 0) {
      stop(func_name, ": none of the input genes were found in the annotation database. ",
           "Check that the gene ID type matches the keytype used in buildAnnot().", call. = FALSE)
    }
  }
  return(input)
}

#' Validate enrichment parameters
#'
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param minSize minimum gene set size
#' @param maxSize maximum gene set size
#' @param minGSSize minimum gene set size for annotation
#' @param maxGSSize maximum gene set size for annotation
#' @param func_name name of calling function
#' @keywords internal
.validateParams <- function(pvalue = 0.05, padj = NULL, minSize = 2, maxSize = 500,
                            minGSSize = 10, maxGSSize = 500, func_name = "enrich") {
  if (!is.null(pvalue) && (pvalue <= 0 || pvalue > 1)) {
    stop(func_name, ": pvalue must be between 0 and 1 (got ", pvalue, ").", call. = FALSE)
  }
  if (!is.null(padj) && (padj <= 0 || padj > 1)) {
    stop(func_name, ": padj must be between 0 and 1 (got ", padj, ").", call. = FALSE)
  }
  if (minSize > maxSize) {
    stop(func_name, ": minSize (", minSize, ") cannot be greater than maxSize (", maxSize, ").", call. = FALSE)
  }
  if (minGSSize > maxGSSize) {
    stop(func_name, ": minGSSize (", minGSSize, ") cannot be greater than maxGSSize (", maxGSSize, ").", call. = FALSE)
  }
  invisible(TRUE)
}

#' Handle empty enrichment results gracefully
#'
#' @param resultFis result data.frame (possibly zero-row)
#' @param pvalue p-value cutoff used
#' @param padj p-adjust cutoff used
#' @param func_name name of calling function
#' @return the input data.frame (may be zero-row)
#' @keywords internal
.handleEmptyResult <- function(resultFis, pvalue = 0.05, padj = NULL, func_name = "enrich") {
  if (nrow(resultFis) == 0) {
    message(func_name, ": no significant terms found. Consider:\n",
            "  - Increasing the pvalue cutoff (current: ",
            if (!is.null(padj)) paste0("padj < ", padj) else paste0("pvalue < ", pvalue), ")\n",
            "  - Adjusting minGSSize/maxGSSize filters\n",
            "  - Checking that input genes overlap with the annotation database")
  }
  return(resultFis)
}
