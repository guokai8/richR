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
#' @details Batch comparison combines each group's already-filtered result rows without
#' applying a second significance cutoff. Successful zero-row results remain
#' in \code{results} and are represented by their group level in
#' \code{combined}; their status is \code{empty}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{results}{Named list of richResult objects}
#'   \item{combined}{Combined data.frame with a \code{group} column,
#'     ready for \code{richCompareDot()}. Its group is a factor whose levels
#'     are the names of all successful results, including empty results.}
#'   \item{status}{One row per input group with character columns \code{group}
#'     and \code{status} (success, empty, failed), integer \code{n_terms}, and
#'     character \code{error}. Failed groups have NA term counts and the
#'     original error message; successful and empty groups have NA errors.}
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
  .validate_named_list(gene_lists, "gene_lists")
  fun <- match.arg(fun)
  n <- length(gene_lists)
  group_names <- names(gene_lists)
  results <- vector("list", n)
  names(results) <- group_names
  statuses <- rep("failed", n)
  n_terms <- rep(NA_integer_, n)
  errors <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    name <- group_names[i]
    if (verbose) message("batchEnrich: [", i, "/", n, "] Processing '", name, "'...")
    tryCatch({
      result <- NULL
      if (fun == "richGO") {
        result <- richGO(gene_lists[[name]], godata = annot, ...)
      } else if (fun == "richKEGG") {
        result <- richKEGG(gene_lists[[name]], kodata = annot, ...)
      } else {
        result <- enrich(gene_lists[[name]], object = annot, ...)
      }
      results[[i]] <- result
      n_terms[[i]] <- as.integer(nrow(result@result))
      statuses[[i]] <- if (n_terms[[i]] > 0L) "success" else "empty"
    }, error = function(e) {
      errors[[i]] <<- conditionMessage(e)
      results[i] <<- list(NULL)
      warning("batchEnrich: failed for '", name, "': ", errors[[i]], call. = FALSE)
    })
  }
  names(results) <- group_names
  results <- results[statuses != "failed"]
  if (length(results) == 0) {
    stop("batchEnrich: all enrichment analyses failed.", call. = FALSE)
  }
  if (verbose) {
    message("batchEnrich: completed ", length(results), "/", n, " analyses successfully.")
  }
  combined_frames <- lapply(names(results), function(name) {
    result <- results[[name]]@result
    result$group <- rep(name, nrow(result))
    result
  })
  combined <- dplyr::bind_rows(combined_frames)
  combined$group <- factor(as.character(combined$group), levels = names(results))
  status <- data.frame(
    group = as.character(group_names), status = as.character(statuses),
    n_terms = as.integer(n_terms), error = as.character(errors),
    stringsAsFactors = FALSE, row.names = NULL
  )
  return(list(results = results, combined = combined, status = status))
}
