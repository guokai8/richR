#' Per-identifier audit of a conversion or an enrichment analysis
#'
#' \code{mappingAudit()} lists what happened to every identifier on its way to
#' the enrichment engine, read only from what was recorded: the mapping record
#' of a \code{\link{mapGenes}()} conversion and the analysis metadata of a
#' result (\code{\link{enrichmentInfo}()}). Nothing is reconstructed from
#' \code{@result}, \code{@gene} or \code{@detail}.
#'
#' @section Availability in existing analyses:
#' Existing enrichment functions do not propagate the record attached by
#' \code{mapGenes()}. Their metadata can still provide a plain audit of the
#' identifiers and background facts they recorded. Joined layouts below are
#' available only when a supported stored record already contains matching
#' conversion evidence; this function does not manufacture that evidence.
#'
#' @section Converted vector:
#' For a vector returned by \code{mapGenes()} the result is its decision table:
#' \code{position}, \code{input_id}, \code{candidates} (a list column),
#' \code{target}, \code{decision}, \code{reason} and, for a ranking,
#' \code{score}; one row per input entry in input order.
#'
#' @section ORA results (richResult):
#' Columns \code{role}, \code{input_id}, \code{mapping}, \code{id},
#' \code{occurrences}, \code{decision}, \code{reason}, in one of two layouts.
#' \describe{
#'   \item{plain layout}{(no mapping record) \code{role = "query"} rows, one per
#'     unique valid query identifier (\code{input_id} = \code{id},
#'     \code{mapping} \code{NA}, \code{occurrences} its recorded count), plus one
#'     \code{invalid} row per missing or blank entry (\code{input_id} and
#'     \code{id} \code{NA}, \code{occurrences} 1).}
#'   \item{joined layout}{(the input carried a matching \code{mapGenes()}
#'     record) one \code{query} row per mapping entry in input order:
#'     \code{input_id} the source identifier, \code{mapping} its mapping
#'     decision, \code{id} its target (\code{NA} when dropped), \code{decision}
#'     \code{not_submitted} for dropped entries and otherwise the analysis
#'     decision of the target, and \code{occurrences} the target's recorded count
#'     on its \code{mapped} row (\code{NA} on \code{merged} and
#'     \code{not_submitted} rows).}
#' }
#' Query decisions: \code{effective} (counted in \code{n}),
#' \code{not_annotated} (no membership in the annotation the analysis tested),
#' \code{outside_background} (not in the supplied universe), \code{invalid}.
#' A stored record explicitly naming \code{universe_policy = "measured"}
#' treats an unannotated query gene in the universe as \code{effective} and
#' one outside it as \code{outside_background}. Reading that record does not
#' add this analysis option to an existing enrichment function.
#' Both layouts add \code{role = "background"} rows, one per unique valid
#' background identifier, \code{effective} (counted in \code{N}) or
#' \code{not_annotated}. The counts reconcile with the engine: the distinct
#' \code{id}s of \code{effective} query rows number \code{testing$effective_n},
#' the \code{effective} background rows number \code{testing$effective_N}, and
#' the non-missing query \code{occurrences} sum to \code{input$original_count}.
#'
#' @section GSEA results (GSEAResult):
#' The same columns plus \code{score} and \code{tied}, all rows
#' \code{role = "ranked"}: one per ranked gene in recorded order (plain layout,
#' \code{occurrences} 1) or one per mapping entry (joined layout,
#' \code{occurrences} 1 on the \code{mapped} row, \code{NA} otherwise).
#' \code{decision} is \code{ranked_annotated}, \code{ranked_not_annotated} (the
#' gene stays in the ranking but belongs to no gene set) or
#' \code{not_submitted}; \code{score} is the gene's value in the result's
#' \code{input} slot (the source score for a \code{not_submitted} row) and
#' \code{tied} is \code{TRUE} when another ranked gene has the same score. The
#' distinct \code{ranked_annotated} ids number \code{input$mapped_count}.
#'
#' @section Older analysis metadata:
#' Metadata without occurrence counts gives the plain layout with
#' \code{occurrences} \code{NA} on identifier rows; metadata without
#' \code{options$universe_policy} is read as \code{"annotated"}. An object
#' without metadata, with the unavailable record, or a vector without a mapping
#' record (or whose \code{richR_mapping} attribute is not a \code{mapGenes()}
#' record) gives \code{NULL} with a warning stating why. A result whose mapping
#' record did not match its input gives the plain layout with a warning.
#'
#' @param x a vector returned by \code{mapGenes()}, a \code{richResult} or a
#'   \code{GSEAResult}.
#' @return a data frame, or \code{NULL} with a warning when no audit is
#'   available.
#' @seealso \code{\link{mapGenes}()}, \code{\link{enrichmentInfo}()}
#' @examples
#' annot <- data.frame(GeneID = paste0("G", 1:20), Term = rep(paste0("T", 1:4), each = 5))
#' res <- enrich(c("G1", "G2", "G3", "G1", "X9"), annot, pvalue = 1, padj = NULL,
#'               minSize = 1, minGSSize = 1)
#' mappingAudit(res)
#' @export
mappingAudit <- function(x) {
  is_rich <- methods::is(x, "richResult")
  is_gsea <- methods::is(x, "GSEAResult")
  if (!is_rich && !is_gsea) {
    record <- if (is.atomic(x)) attr(x, "richR_mapping", exact = TRUE)
    if (is.null(record)) return(.audit_unavailable("the object carries no mapping record"))
    if (!.mapping_record_ok(record)) return(.audit_unavailable("the mapping record is malformed"))
    return(.mapping_decision_frame(record))
  }
  m <- attr(x, "richR_metadata", exact = TRUE)
  if (!is.list(m)) return(.audit_unavailable("metadata not recorded"))
  if (!isTRUE(m$available)) {
    return(.audit_unavailable(m$reason %||% "metadata unavailable"))
  }
  mapping <- m$input$mapping
  if (is.list(mapping) && identical(mapping$available, FALSE)) {
    warning("mappingAudit(): the mapping record is unavailable (",
            mapping$reason %||% "no reason recorded",
            "); showing the analysis input only.", call. = FALSE)
    mapping <- NULL
  }
  if (is_gsea) .audit_gsea(x, m, mapping) else .audit_ora(m, mapping)
}

## The record's decision table (a plain column list) as a data frame with the
## same columns in the same order, candidates as a list column, automatic row
## names.
.mapping_decision_frame <- function(record) {
  structure(record$decisions, names = names(record$decisions),
            row.names = .set_row_names(length(record$decisions[[1L]])),
            class = "data.frame")
}

## TRUE when record carries a decision table of equal-length columns with the
## character decision and target columns the audit reads.
.mapping_record_ok <- function(record) {
  d <- if (is.list(record)) record$decisions
  is.list(d) && length(d) > 0L && is.character(d$decision) && is.character(d$target) &&
    length(unique(lengths(d))) == 1L
}

.audit_unavailable <- function(why) {
  warning("mappingAudit(): mapping audit unavailable: ", why, ".", call. = FALSE)
  NULL
}

## The record of x when its kept targets (and scores) equal x; an unavailable
## marker, with a warning, when they do not; NULL when x carries none.
.mapping_for_input <- function(x) {
  record <- if (is.atomic(x)) attr(x, "richR_mapping", exact = TRUE)
  if (is.null(record)) return(NULL)
  d <- if (is.list(record)) record$decisions
  ok <- is.list(d) && is.character(d$decision) && is.character(d$target)
  if (ok) {
    keep <- d$decision %in% "mapped"
    ok <- if (!is.null(d$score)) {
      is.numeric(x) && identical(names(x), d$target[keep]) &&
        identical(unname(as.double(x)), d$score[keep])
    } else {
      is.character(x) && identical(unname(as.character(x)), d$target[keep])
    }
  }
  if (ok) return(record)
  reason <- "mapping record does not match the submitted input"
  warning(reason, call. = FALSE)
  list(available = FALSE, reason = reason)
}

## Analysis decision of each recorded query identifier, read from metadata.
.ora_query_decisions <- function(m) {
  inp <- m$input
  ids <- as.character(inp$ids)
  policy <- m$options$universe_policy %||% "annotated"
  explicit <- identical(m$background$mode, "explicit")
  annotated <- ids %in% inp$mapped_ids
  effective <- ids %in% inp$effective_ids
  in_universe_unannotated <- ids %in% m$background$unmapped_ids
  decision <- character(length(ids))
  reason <- character(length(ids))
  decision[effective & annotated] <- "effective"
  reason[effective & annotated] <- "annotated and in the background; counted in n"
  e_un <- effective & !annotated
  decision[e_un] <- "effective"
  reason[e_un] <- paste("not annotated; in the supplied universe and counted in n",
                        "under universe_policy = \"measured\"")
  out_ann <- !effective & annotated
  decision[out_ann] <- "outside_background"
  reason[out_ann] <- "annotated but not in the supplied universe"
  rest <- !effective & !annotated
  if (policy == "measured") {
    decision[rest] <- "outside_background"
    reason[rest] <- "not annotated and not in the supplied universe"
  } else {
    decision[rest] <- "not_annotated"
    reason[rest] <- if (!explicit) {
      "not annotated in the tested annotation"
    } else {
      ifelse(in_universe_unannotated[rest],
             "not annotated in the tested annotation; in the supplied universe",
             "not annotated in the tested annotation; outside the supplied universe")
    }
  }
  occ <- inp$occurrences
  occ <- if (is.null(occ)) rep(NA_integer_, length(ids)) else as.integer(occ)
  list(ids = ids, decision = decision, reason = reason, occurrences = occ)
}

.audit_frame <- function(role, input_id, mapping, id, occurrences, decision, reason) {
  data.frame(role = role, input_id = input_id, mapping = mapping, id = id,
             occurrences = as.integer(occurrences), decision = decision,
             reason = reason, stringsAsFactors = FALSE)
}

.audit_ora <- function(m, mapping) {
  qd <- .ora_query_decisions(m)
  if (is.null(mapping)) {
    n_invalid <- as.integer(m$input$invalid_count %||% 0L)
    if (is.na(n_invalid)) n_invalid <- 0L
    query <- .audit_frame(
      "query", c(qd$ids, rep(NA_character_, n_invalid)), NA_character_,
      c(qd$ids, rep(NA_character_, n_invalid)),
      c(qd$occurrences, rep(1L, n_invalid)),
      c(qd$decision, rep("invalid", n_invalid)),
      c(qd$reason, rep("missing or blank identifier", n_invalid)))
  } else {
    d <- mapping$decisions
    at <- match(d$target, qd$ids)
    kept <- !is.na(d$target)
    decision <- ifelse(kept, qd$decision[at], "not_submitted")
    reason <- ifelse(kept, qd$reason[at], paste0("not submitted (", d$decision, "): ", d$reason))
    occ <- ifelse(d$decision == "mapped", qd$occurrences[at], NA_integer_)
    query <- .audit_frame("query", d$input_id, d$decision, d$target, occ, decision, reason)
  }
  bg <- m$background
  b_ids <- unique(c(as.character(bg$effective_ids), as.character(bg$unmapped_ids)))
  b_eff <- b_ids %in% bg$effective_ids
  b_unann <- b_ids %in% bg$unmapped_ids
  b_reason <- ifelse(b_eff & !b_unann, "annotated; counted in N",
              ifelse(b_eff, paste("not annotated; counted in N under",
                                  "universe_policy = \"measured\""),
                     "not annotated in the tested annotation; not counted in N"))
  background <- .audit_frame("background", b_ids, NA_character_, b_ids, NA_integer_,
                             ifelse(b_eff, "effective", "not_annotated"), b_reason)
  out <- rbind(query, background)
  rownames(out) <- NULL
  out
}

.audit_gsea <- function(x, m, mapping) {
  inp <- m$input
  ids <- as.character(inp$ids)
  ranked <- x@input
  tie <- duplicated(ranked) | duplicated(ranked, fromLast = TRUE)
  annotated <- ids %in% inp$mapped_ids
  decision <- ifelse(annotated, "ranked_annotated", "ranked_not_annotated")
  reason <- ifelse(annotated, "in at least one gene set",
                   "in no gene set; stays in the ranking")
  pos <- match(ids, names(ranked))
  score <- unname(as.double(ranked[pos]))
  tied <- unname(tie[pos])
  occ <- if (is.null(inp$occurrences)) rep(NA_integer_, length(ids)) else as.integer(inp$occurrences)
  if (is.null(mapping)) {
    out <- .audit_frame("ranked", ids, NA_character_, ids, occ, decision, reason)
    out$score <- score
    out$tied <- tied
    return(out)
  }
  d <- mapping$decisions
  at <- match(d$target, ids)
  kept <- !is.na(d$target)
  out <- .audit_frame(
    "ranked", d$input_id, d$decision, d$target,
    ifelse(d$decision == "mapped", occ[at], NA_integer_),
    ifelse(kept, decision[at], "not_submitted"),
    ifelse(kept, reason[at], paste0("not submitted (", d$decision, "): ", d$reason)))
  out$score <- ifelse(kept, score[at], as.double(d$score))
  out$tied <- ifelse(kept, tied[at], NA)
  out
}
