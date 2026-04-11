## ----------------------------------------------------------------
##  show methods
## ----------------------------------------------------------------

#' Show methods for richR S4 classes
#'
#' Display a concise summary when an object is printed at the console.
#'
#' @param object An S4 object of class richResult, GSEAResult, or Annot
#' @name show-methods
#' @aliases show,richResult-method
#' @rdname show-methods
setMethod("show", signature(object = "richResult"), function(object) {
  n_sig <- nrow(object@result)
  cat("richResult object\n")
  cat("  Organism:    ", if (length(object@organism) && nchar(object@organism)) object@organism else "(not set)", "\n")
  cat("  Ontology:    ", if (length(object@ontology) && nchar(object@ontology)) object@ontology else "(not set)", "\n")
  cat("  Gene count:  ", object@genenumber, "input genes\n")
  cat("  Significant: ", n_sig, "terms")
  if (length(object@pvalueCutoff) && object@pvalueCutoff > 0) {
    cat(" (pvalue <", object@pvalueCutoff, ")")
  }
  cat("\n")
  if (n_sig > 0) {
    top_n <- min(5, n_sig)
    top <- object@result[1:top_n, c("Term", "Pvalue", "Padj", "Significant", "Annotated"), drop = FALSE]
    top$Pvalue <- formatC(top$Pvalue, format = "e", digits = 2)
    top$Padj <- formatC(as.numeric(top$Padj), format = "e", digits = 2)
    cat("\n  Top", top_n, "terms:\n")
    print(top, row.names = FALSE, right = FALSE)
    if (n_sig > top_n) cat("  ... with", n_sig - top_n, "more terms\n")
  }
  invisible(object)
})

#' @rdname show-methods
setMethod("show", signature(object = "GSEAResult"), function(object) {
  n_sig <- nrow(object@result)
  cat("GSEAResult object\n")
  cat("  Organism:    ", if (length(object@organism) && nchar(object@organism)) object@organism else "(not set)", "\n")
  cat("  Ontology:    ", if (length(object@ontology) && nchar(object@ontology)) object@ontology else "(not set)", "\n")
  cat("  Gene count:  ", object@genenumber, "ranked genes\n")
  cat("  Significant: ", n_sig, "pathways\n")
  if (n_sig > 0) {
    top_n <- min(5, n_sig)
    cols <- intersect(c("pathway", "pval", "padj", "NES", "size"), colnames(object@result))
    top <- object@result[1:top_n, cols, drop = FALSE]
    if ("pval" %in% cols) top$pval <- formatC(top$pval, format = "e", digits = 2)
    if ("padj" %in% cols) top$padj <- formatC(top$padj, format = "e", digits = 2)
    if ("NES" %in% cols) top$NES <- round(top$NES, 3)
    cat("\n  Top", top_n, "pathways:\n")
    print(top, row.names = FALSE, right = FALSE)
    if (n_sig > top_n) cat("  ... with", n_sig - top_n, "more pathways\n")
  }
  invisible(object)
})

#' @rdname show-methods
#' @aliases show,Annot-method
setMethod("show", signature(object = "Annot"), function(object) {
  n_genes <- length(unique(object@annot[, 1]))
  n_terms <- length(unique(object@annot[, 2]))
  cat("Annot object\n")
  cat("  Species: ", if (nchar(object@species)) object@species else "(not set)", "\n")
  cat("  Type:    ", if (nchar(object@anntype)) object@anntype else "(not set)", "\n")
  cat("  Keytype: ", if (nchar(object@keytype)) object@keytype else "(not set)", "\n")
  cat("  Genes:   ", n_genes, "\n")
  cat("  Terms:   ", n_terms, "\n")
  invisible(object)
})

## ----------------------------------------------------------------
##  getGenes
## ----------------------------------------------------------------

#' Extract genes from enrichment results
#'
#' Returns genes associated with a specific term, or all input genes if no term
#' is specified.
#'
#' @param x A richResult or GSEAResult object
#' @param term Optional term ID (Annot column) or term name to extract genes for.
#'   If NULL, returns all input genes from the analysis.
#' @param sep Separator used in GeneID column (default: ",")
#' @return Character vector of gene IDs
#'
#' @examples
#' \dontrun{
#' res <- richGO(genes, godata = hsago, ontology = "BP")
#' # All input genes
#' getGenes(res)
#' # Genes in a specific term
#' getGenes(res, "GO:0006915")
#' }
#' @export
#' @author Junguk Hur
getGenes <- function(x, term = NULL, sep = ",") {
  if (inherits(x, "richResult") || inherits(x, "GSEAResult")) {
    if (is.null(term)) {
      return(x@gene)
    }
    res_df <- x@result
    if (inherits(x, "richResult")) {
      # Match by Annot ID or Term name
      idx <- which(res_df$Annot == term | res_df$Term == term)
      if (length(idx) == 0) {
        stop("Term '", term, "' not found in results.", call. = FALSE)
      }
      gene_str <- as.character(res_df$GeneID[idx[1]])
      return(unique(trimws(strsplit(gene_str, sep)[[1]])))
    }
    if (inherits(x, "GSEAResult")) {
      idx <- which(res_df$pathway == term)
      if (length(idx) == 0) {
        stop("Pathway '", term, "' not found in GSEA results.", call. = FALSE)
      }
      le <- res_df$leadingEdge[idx[1]]
      if (is.list(le)) return(unique(le[[1]]))
      return(unique(strsplit(as.character(le), sep)[[1]]))
    }
  }
  stop("x must be a richResult or GSEAResult object.", call. = FALSE)
}

## ----------------------------------------------------------------
##  as.data.frame methods
## ----------------------------------------------------------------

##' @method as.data.frame Annot
##' @export
as.data.frame.Annot <- function(x, ...) {
  as.data.frame(x@annot)
}

##' @method as.data.frame richResult
##' @export
as.data.frame.richResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}

##' @method as.data.frame GSEAResult
##' @export
as.data.frame.GSEAResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}

## ----------------------------------------------------------------
##  setAs coercion methods
## ----------------------------------------------------------------

setAs(from = "richResult", to = "data.frame", def = function(from) {
  as.data.frame(from@result)
})

setAs(from = "GSEAResult", to = "data.frame", def = function(from) {
  as.data.frame(from@result)
})

setAs(from = "Annot", to = "data.frame", def = function(from) {
  as.data.frame(from@annot)
})

setAs(from = "data.frame", to = "Annot", def = function(from) {
  keytype <- character()
  species <- character()
  anntype <- character()
  GeneID <- as.vector(from[, 1])
  Term <- as.vector(from[, 2])
  Annot <- from$Annot
  annot <- data.frame(GeneID, Term, Annot)
  new("Annot",
      species = species,
      anntype = anntype,
      keytype = keytype,
      annot = annot
  )
})

## ----------------------------------------------------------------
##  row.names / rownames / names / colnames
## ----------------------------------------------------------------

#' S4 row- and column-related methods for Annot, richResult, GSEAResult
#' @param x An object of class Annot, richResult, or GSEAResult
#' @name s4-accessors
#' @rdname s4-accessors
NULL

#' @rdname s4-accessors
#' @aliases row.names,Annot-method
setMethod("row.names", signature(x = "Annot"), function(x) row.names(x@annot))

#' @rdname s4-accessors
#' @aliases row.names,richResult-method
setMethod("row.names", signature(x = "richResult"), function(x) row.names(x@result))

#' @rdname s4-accessors
#' @aliases row.names,GSEAResult-method
setMethod("row.names", signature(x = "GSEAResult"), function(x) row.names(x@result))

#' @rdname s4-accessors
#' @aliases rownames,Annot-method
setMethod("rownames", signature(x = "Annot"), function(x) rownames(x@annot))

#' @rdname s4-accessors
#' @aliases rownames,richResult-method
setMethod("rownames", signature(x = "richResult"), function(x) rownames(x@result))

#' @rdname s4-accessors
#' @aliases rownames,GSEAResult-method
setMethod("rownames", signature(x = "GSEAResult"), function(x) rownames(x@result))

#' @rdname s4-accessors
#' @aliases names,Annot-method
setMethod("names", signature(x = "Annot"), function(x) names(x@annot))

#' @rdname s4-accessors
#' @aliases names,richResult-method
setMethod("names", signature(x = "richResult"), function(x) names(x@result))

#' @rdname s4-accessors
#' @aliases names,GSEAResult-method
setMethod("names", signature(x = "GSEAResult"), function(x) names(x@result))

#' @rdname s4-accessors
#' @aliases colnames,Annot-method
setMethod("colnames", signature(x = "Annot"), function(x) colnames(x@annot))

#' @rdname s4-accessors
#' @aliases colnames,richResult-method
setMethod("colnames", signature(x = "richResult"), function(x) colnames(x@result))

#' @rdname s4-accessors
#' @aliases colnames,GSEAResult-method
setMethod("colnames", signature(x = "GSEAResult"), function(x) colnames(x@result))

## ----------------------------------------------------------------
##  head / tail / dim
## ----------------------------------------------------------------

#' S4 head, tail, dim methods
#' @param x An S4 object
#' @param n Number of rows
#' @param ... Further arguments
#' @name utilities-methods
#' @rdname utilities-methods
NULL

#' @rdname utilities-methods
#' @aliases head,Annot-method
setMethod("head", signature(x = "Annot"), function(x, n = 6L, ...) {
  cat("=== species is:", x@species, "and Annotation is", x@anntype,
      "keytype is", x@keytype, "===\n")
  utils::head(x@annot, n, ...)
})

#' @rdname utilities-methods
#' @aliases head,richResult-method
setMethod("head", signature(x = "richResult"), function(x, n = 6L, ...) {
  cat("=== Total significant terms is:", dim(x@result), "===\n")
  utils::head(x@result, n, ...)
})

#' @rdname utilities-methods
#' @aliases head,GSEAResult-method
setMethod("head", signature(x = "GSEAResult"), function(x, n = 6L, ...) {
  cat("=== Total significant terms is:", dim(x@result), "===\n")
  utils::head(x@result, n, ...)
})

#' @rdname utilities-methods
#' @aliases tail,Annot-method
setMethod("tail", signature(x = "Annot"), function(x, n = 6L, ...) {
  cat("=== species is:", x@species, "and Annotation is", x@anntype,
      "keytype is", x@keytype, "===\n")
  utils::tail(x@annot, n, ...)
})

#' @rdname utilities-methods
#' @aliases tail,richResult-method
setMethod("tail", signature(x = "richResult"), function(x, n = 6L, ...) {
  cat("=== Total significant terms is:", dim(x@result), "===\n")
  utils::tail(x@result, n, ...)
})

#' @rdname utilities-methods
#' @aliases tail,GSEAResult-method
setMethod("tail", signature(x = "GSEAResult"), function(x, n = 6L, ...) {
  cat("=== Total significant terms is:", dim(x@result), "===\n")
  utils::tail(x@result, n, ...)
})

#' @rdname utilities-methods
#' @aliases dim,Annot-method
setMethod("dim", signature(x = "Annot"), function(x) dim(x@annot))

#' @rdname utilities-methods
#' @aliases dim,richResult-method
setMethod("dim", signature(x = "richResult"), function(x) dim(x@result))

#' @rdname utilities-methods
#' @aliases dim,GSEAResult-method
setMethod("dim", signature(x = "GSEAResult"), function(x) dim(x@result))

## ----------------------------------------------------------------
##  summary
## ----------------------------------------------------------------

#' S4 summary method for richResult objects
#' @param object richResult object
#' @param ... Not used
#' @return Prints summary to console
#' @export
#' @aliases summary,richResult-method
setMethod("summary", signature = "richResult", definition = function(object, ...) {
  cat("Total input genes is:", length(object@gene),
      "and significant biological term is:", nrow(object@result), "\n")
})

#' S4 summary method for GSEAResult objects
#' @param object GSEAResult object
#' @param ... Not used
#' @return Prints summary to console
#' @export
#' @aliases summary,GSEAResult-method
setMethod("summary", signature = "GSEAResult", definition = function(object, ...) {
  n_sig <- sum(object@result$padj < 0.05, na.rm = TRUE)
  cat("Total significant biological term is:", n_sig, "\n")
})

## ----------------------------------------------------------------
##  [ subsetting
## ----------------------------------------------------------------

#' S4 subsetting methods
#' @param x An S4 object
#' @param i,j indices
#' @param ... further arguments
#' @param drop logical
#' @name subset-methods
#' @rdname subset-methods
NULL

#' @rdname subset-methods
#' @aliases [,Annot-method
setMethod("[", signature(x = "Annot", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) x@annot[i, j, drop = drop])

#' @rdname subset-methods
#' @aliases [,richResult-method
setMethod("[", signature(x = "richResult", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) x@result[i, j, drop = drop])

#' @rdname subset-methods
#' @aliases [,GSEAResult-method
setMethod("[", signature(x = "GSEAResult", i = "ANY", j = "ANY", drop = "ANY"),
          function(x, i, j, ..., drop = TRUE) x@result[i, j, drop = drop])

## ----------------------------------------------------------------
##  $ accessor
## ----------------------------------------------------------------

#' S4 dollar accessor methods
#' @param x An S4 object
#' @param name column name
#' @name dollar-methods
#' @rdname dollar-methods
NULL

#' @rdname dollar-methods
#' @aliases $,Annot-method
setMethod("$", signature(x = "Annot"), function(x, name) x@annot[[name]])

#' @rdname dollar-methods
#' @aliases $,richResult-method
setMethod("$", signature(x = "richResult"), function(x, name) x@result[[name]])

#' @rdname dollar-methods
#' @aliases $,GSEAResult-method
setMethod("$", signature(x = "GSEAResult"), function(x, name) x@result[[name]])

## ----------------------------------------------------------------
##  result() and detail() S4 methods
## ----------------------------------------------------------------

#' S4 methods for result() and detail()
#' @param x richResult or GSEAResult object
#' @return data.frame
#' @name result_detail_methods
#' @rdname result_detail_methods
NULL

#' @rdname result_detail_methods
#' @aliases result,richResult-method
#' @export
setMethod("result", signature = "richResult", function(x) as.data.frame(x@result))

#' @rdname result_detail_methods
#' @aliases result,GSEAResult-method
#' @export
setMethod("result", signature = "GSEAResult", function(x) as.data.frame(x@result))

#' @rdname result_detail_methods
#' @aliases detail,richResult-method
#' @export
setMethod("detail", signature = "richResult", function(x) as.data.frame(x@detail))

#' @rdname result_detail_methods
#' @aliases detail,GSEAResult-method
#' @export
setMethod("detail", signature = "GSEAResult", function(x) {
  message("GSEAResult does not have a detail slot; returning result instead")
  as.data.frame(x@result)
})

## ----------------------------------------------------------------
##  rbind methods
## ----------------------------------------------------------------

#' rbind for richResult objects
#' @importFrom S4Vectors bindROWS
#' @export
#' @author Kai Guo
rbind.richResult <- function(...) {
  objects <- list(...)
  objects <- lapply(objects, as.data.frame)
  bindROWS(objects[[1L]], objects = objects[-1L])
}

#' rbind for GSEAResult objects
#' @importFrom S4Vectors bindROWS
#' @export
#' @author Kai Guo
rbind.GSEAResult <- function(...) {
  objects <- list(...)
  objects <- lapply(objects, as.data.frame)
  bindROWS(objects[[1L]], objects = objects[-1L])
}
