#' Tested terms of an enrichment result
#'
#' Returns the tested table of an enrichment result: every term the analysis
#' tested, not only the terms it reported in \code{result(x)}.
#'
#' For an ORA result (\code{richResult}) the table holds one row for every
#' term of the tested family whose overlap with the query is greater than
#' zero. This upstream engine adjusts positive-overlap terms before size
#' reporting filters; its recorded zero-overlap count is zero in
#' \code{enrichmentInfo(x)$testing$zero_overlap_count} but not listed. For a
#' GSEA result (\code{GSEAResult}) it holds one row for every pathway fgsea
#' returned. The columns are those of \code{x@result}, with the values the
#' analysis computed, plus a logical column \code{reported} that is
#' \code{TRUE} exactly for the terms the analysis reported. Terms removed by
#' the size filters, \code{minSize}, \code{keepRich} or the significance
#' cutoffs stay in the table with \code{reported = FALSE}.
#'
#' Objects created by a richR version without the tested table, by \code{as(<data.frame>, "richResult")}, by \code{richDAVID()} or
#' by \code{new()} carry no recorded tested table. For them
#' \code{testedTerms()} returns \code{NULL} with a warning; the table is never
#' reconstructed from \code{result(x)}.
#'
#' @param x a \code{richResult} or \code{GSEAResult} object
#' @return a data.frame (zero rows when the analysis tested nothing), or
#'   \code{NULL} with a warning when the object carries no recorded tested
#'   table
#' @seealso \code{\link{enrichmentInfo}}
#' @export
testedTerms <- function(x) {
  if (!(methods::is(x, "richResult") || methods::is(x, "GSEAResult"))) {
    stop("testedTerms() supports richResult or GSEAResult objects.",
         call. = FALSE)
  }
  reason <- .tested_unavailable_reason(x)
  if (!is.null(reason)) {
    warning("tested table unavailable: ", reason, call. = FALSE)
    return(NULL)
  }
  x@tested
}

#' Why a result object carries no recorded tested table
#'
#' The single availability rule for the \code{tested} slot: every reader of
#' the slot goes through here, so an object deserialized from a richR version
#' without the slot never reaches \code{x@tested}. A recorded tested table
#' always has columns (at least the identifier and \code{reported}); the slot
#' prototype, a zero-column data.frame, means "not recorded".
#' @param x a richResult or GSEAResult object
#' @return NULL when the tested table is recorded, otherwise the reason
#' @keywords internal
#' @noRd
.tested_unavailable_reason <- function(x) {
  if (!methods::.hasSlot(x, "tested")) {
    return(paste("the object has no tested slot (created by a richR version",
                 "without tested-evidence support); it is not reconstructed from @result"))
  }
  if (ncol(x@tested) == 0L) {
    ## Neutral on provenance: the object may come from as(), richDAVID() or
    ## new(), or be a legacy result read back by readRich().
    return(paste("no tested table was recorded for this object;",
                 "it is not reconstructed from @result"))
  }
  NULL
}

#' Does a result object carry a recorded tested table?
#' @param x a richResult or GSEAResult object
#' @return TRUE or FALSE
#' @keywords internal
#' @noRd
# @MX:ANCHOR: [AUTO] Shared tested-evidence availability contract.
# @MX:REASON: Binding, storage and testedTerms must agree on historical absence.
.has_tested <- function(x) {
  is.null(.tested_unavailable_reason(x))
}
