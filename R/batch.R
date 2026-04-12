#' Batch enrichment analysis across multiple gene lists
#'
#' Runs enrichment analysis on a named list of gene vectors and returns
#' combined results ready for \code{richCompareDot()} visualization.
#'
#' @param gene_lists Named list of character vectors (gene IDs per group)
#' @param annot Annot object or annotation data.frame
#' @param fun Enrichment function to use: "richGO", "richKEGG", or "enrich"
#'   (default: "enrich")
#' @param verbose Logical, print progress messages (default: TRUE)
#' @param ... Additional arguments passed to the enrichment function
#'   (e.g., ontology, pvalue, padj, minSize, maxSize)
#'
#' @return A list with two elements:
#' \describe{
#'   \item{results}{Named list of richResult objects}
#'   \item{combined}{Combined data.frame with a \code{group} column,
#'     ready for \code{richCompareDot()}}
#' }
#'
#' @examples
#' \dontrun{
#' hsako <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")
#' gene_lists <- list(
#'   "Treatment" = sample(unique(hsako$GeneID), 500),
#'   "Control"   = sample(unique(hsako$GeneID), 500)
#' )
#'
#' # Run batch enrichment
#' batch <- batchEnrich(gene_lists, hsako, fun = "richKEGG")
#'
#' # Access individual results
#' head(result(batch$results[["Treatment"]]))
#'
#' # Visualize comparison
#' richCompareDot(batch$combined)
#' }
#' @export
#' @author Junguk Hur
batchEnrich <- function(gene_lists, annot, fun = c("enrich", "richGO", "richKEGG"),
                        verbose = TRUE, ...) {
  fun <- match.arg(fun)
  if (!is.list(gene_lists) || is.null(names(gene_lists))) {
    stop("gene_lists must be a named list of character vectors.", call. = FALSE)
  }
  n <- length(gene_lists)
  group_names <- names(gene_lists)
  results <- list()
  for (i in seq_len(n)) {
    name <- group_names[i]
    if (verbose) message("batchEnrich: [", i, "/", n, "] Processing '", name, "'...")
    tryCatch({
      if (fun == "richGO") {
        results[[name]] <- richGO(gene_lists[[name]], godata = annot, ...)
      } else if (fun == "richKEGG") {
        results[[name]] <- richKEGG(gene_lists[[name]], kodata = annot, ...)
      } else {
        results[[name]] <- enrich(gene_lists[[name]], object = annot, ...)
      }
    }, error = function(e) {
      warning("batchEnrich: failed for '", name, "': ", e$message, call. = FALSE)
      results[[name]] <<- NULL
    })
  }
  # Remove failed results
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) {
    stop("batchEnrich: all enrichment analyses failed.", call. = FALSE)
  }
  if (verbose) {
    message("batchEnrich: completed ", length(results), "/", n, " analyses successfully.")
  }
  # Create combined data.frame
  combined <- compareResult(results)
  return(list(results = results, combined = combined))
}
