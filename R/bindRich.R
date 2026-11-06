#' Bind enrichment results into one table with an explicit state per term
#'
#' Combines several enrichment results of one class into one long table with
#' one row per group and per identifier, stating for each row whether the
#' group reported the term, tested it without reporting it, did not test it,
#' or cannot tell. Unlike \code{\link{compareResult}}, nothing is filtered and
#' a term absent from a group is never silently dropped or filled with zeros.
#'
#' The rows cover the union of the identifiers (\code{Annot} for ORA,
#' \code{pathway} for GSEA) found in any group's \code{result(x)} or tested
#' table (\code{\link{testedTerms}}). The columns are, in this order:
#' \describe{
#'   \item{\code{group}}{a factor whose levels are the group names in input
#'     order;}
#'   \item{the identifier}{\code{Annot} (ORA) or \code{pathway} (GSEA);}
#'   \item{\code{Term}}{the group's own label for the term;}
#'   \item{\code{state}}{the inferential state described below;}
#'   \item{\code{execution}}{\code{"success"} (the group's \code{result(x)}
#'     has at least one row), \code{"empty"} (it has zero rows) or
#'     \code{"failed"} (a \code{NULL} entry, as \code{batchEnrich()} returns
#'     for a failed group);}
#'   \item{the remaining columns}{the union of the other columns of every
#'     group's result and tested table, in order of first appearance (groups
#'     in input order, the result before the tested table), typed as at their
#'     first appearance; integer is widened to double when groups disagree,
#'     and a cell a group lacks is a typed \code{NA}. The tested table's
#'     \code{reported} column is not carried over.}
#' }
#'
#' \code{state} is \code{NA} on every row of a failed group and otherwise one
#' of:
#' \describe{
#'   \item{\code{"reported"}}{the identifier is in the group's current
#'     \code{result(x)};}
#'   \item{\code{"tested_not_reported"}}{it is in the group's tested table but
#'     not in its current result (filtered by a cutoff, a size rule, or a
#'     dplyr verb applied to the object);}
#'   \item{\code{"not_tested"}}{it is absent from a group whose tested table
#'     is recorded and complete (\code{zero_overlap_count = 0});}
#'   \item{\code{"untested_or_zero_overlap"}}{it is absent from a group whose
#'     tested table is recorded and whose metadata counts zero-overlap tested
#'     terms that the table does not list in a supported stored record: the
#'     term was either not tested or tested with zero overlap;}
#'   \item{\code{"unavailable"}}{it is absent from the result of a group
#'     without a recorded tested table (objects from versions without tested-evidence support,
#'     \code{as()}, \code{richDAVID()}), or of a group whose metadata is
#'     missing or lacks \code{zero_overlap_count}: what was tested cannot be
#'     known and is not reconstructed.}
#' }
#'
#' \code{Term} and the statistics of a row come only from that group's own
#' record: its result for \code{"reported"}, its tested table for
#' \code{"tested_not_reported"}, and \code{NA} otherwise. No cutoff is applied,
#' no p-value is re-adjusted, and no statistic is computed. A grouped
#' \code{result(x)} (after \code{group_by()}) is read as its plain rows.
#'
#' @param x a non-empty list of \code{richResult} objects, or of
#'   \code{GSEAResult} objects, with \code{NULL} for failed groups. Unnamed
#'   lists receive the labels \code{Group 1}, \code{Group 2}, ...; explicit
#'   names must be nonblank and unique, as in \code{\link{compareResult}}. ORA
#'   and GSEA results cannot be mixed: their statistics are not comparable.
#' @return a data.frame as described above
#' @seealso \code{\link{testedTerms}}, \code{\link{compareResult}}
#' @examples
#' \dontrun{
#' res <- batchEnrich(list(S1 = gene1, S2 = gene2), annot = hsako)
#' bound <- bindRich(res$results)
#' table(bound$group, bound$state)
#' }
#' @export
bindRich <- function(x) {
  if (!is.list(x) || is.data.frame(x) || length(x) == 0L) {
    stop("x must be a non-empty list of enrichment results.", call. = FALSE)
  }
  if (is.null(names(x))) names(x) <- paste("Group", seq_along(x))
  .validate_named_list(x, "x")
  groups <- names(x)

  failed <- vapply(x, is.null, logical(1))
  is_ora <- vapply(x, methods::is, logical(1), "richResult")
  is_gsea <- vapply(x, methods::is, logical(1), "GSEAResult")
  bad <- !(failed | is_ora | is_gsea)
  if (any(bad)) {
    found <- vapply(x[bad], function(v) class(v)[[1]], character(1))
    stop("bindRich: x group(s) ", .bind_quote(groups[bad]),
         " must be a richResult, GSEAResult, or NULL (found ",
         paste(found, collapse = ", "), ").", call. = FALSE)
  }
  if (all(failed)) {
    stop("bindRich: every group of x is NULL (failed); there is no result ",
         "to bind.", call. = FALSE)
  }
  if (any(is_ora) && any(is_gsea)) {
    stop("bindRich() cannot bind ORA and GSEA results together: ",
         "ORA (richResult) group(s) ", .bind_quote(groups[is_ora]),
         "; GSEA (GSEAResult) group(s) ", .bind_quote(groups[is_gsea]),
         ". Bind each class separately.", call. = FALSE)
  }
  id <- if (any(is_ora)) "Annot" else "pathway"

  records <- lapply(seq_along(x), function(i) {
    if (failed[[i]]) return(NULL)
    .bind_record(x[[i]], groups[[i]], id)
  })

  ids <- unique(as.character(unlist(lapply(records, function(r) {
    c(as.character(r$result[[id]]), as.character(r$tested[[id]]))
  }), use.names = FALSE)))
  n_id <- length(ids)
  n <- n_id * length(groups)

  # Value columns in order of first appearance; Term leads.
  value_cols <- unique(c("Term", unlist(lapply(records, function(r) {
    c(names(r$result), setdiff(names(r$tested), "reported"))
  }), use.names = FALSE)))
  value_cols <- setdiff(value_cols, id)
  clash <- intersect(value_cols, c("group", "state", "execution"))
  if (length(clash)) {
    stop("bindRich: result column(s) ", .bind_quote(clash),
         " clash with the columns bindRich() adds.", call. = FALSE)
  }

  # Per row: which group, which source row, which state.
  state <- rep(NA_character_, n)
  execution <- character(n)
  res_pos <- tst_pos <- rep(NA_integer_, n)
  for (i in seq_along(groups)) {
    rows <- (i - 1L) * n_id + seq_len(n_id)
    r <- records[[i]]
    if (is.null(r)) {
      execution[rows] <- "failed"
      next
    }
    execution[rows] <- if (nrow(r$result) > 0L) "success" else "empty"
    rp <- match(ids, as.character(r$result[[id]]))
    tp <- if (is.null(r$tested)) {
      rep(NA_integer_, n_id)
    } else {
      match(ids, as.character(r$tested[[id]]))
    }
    tp[!is.na(rp)] <- NA_integer_
    st <- rep(r$absent, n_id)
    st[!is.na(tp)] <- "tested_not_reported"
    st[!is.na(rp)] <- "reported"
    state[rows] <- st
    res_pos[rows] <- rp
    tst_pos[rows] <- tp
  }
  group_of <- rep(seq_along(groups), each = n_id)

  values <- lapply(value_cols, function(col) {
    .bind_column(col, records, groups, group_of, res_pos, tst_pos, n)
  })
  names(values) <- value_cols

  cols <- c(
    list(group = factor(rep(groups, each = n_id), levels = groups)),
    stats::setNames(list(rep(ids, times = length(groups))), id),
    values["Term"],
    list(state = state, execution = execution),
    values[setdiff(value_cols, "Term")]
  )
  structure(cols, names = names(cols), row.names = .set_row_names(n),
            class = "data.frame")
}

#' One group's record for bindRich(): plain result rows, the tested table
#' when recorded, and the state of an identifier absent from both
#' @param obj a richResult or GSEAResult object
#' @param group the group name, for messages
#' @param id the identifier column
#' @return a list with result, tested (or NULL) and absent
#' @keywords internal
#' @noRd
.bind_record <- function(obj, group, id) {
  result <- as.data.frame(obj@result)
  if (identical(id, "pathway") && "pathway" %in% names(result) && !"Term" %in% names(result)) {
    result$Term <- as.character(result$pathway)
  }
  has_tested <- .has_tested(obj)
  tested <- if (has_tested) as.data.frame(obj@tested) else NULL
  if (!is.null(tested) && identical(id, "pathway") && "pathway" %in% names(tested) && !"Term" %in% names(tested)) {
    tested$Term <- as.character(tested$pathway)
  }
  tables <- list(`@result` = result, `tested table` = tested)
  for (slot in names(tables)) {
    tab <- tables[[slot]]
    if (is.null(tab)) next
    missing <- setdiff(c(id, "Term"), names(tab))
    if (length(missing)) {
      stop("bindRich: x group '", group, "' ", slot, " lacks column(s) ",
           .bind_quote(missing), "; bindRich() binds by identifier and ",
           "label and cannot bind it.", call. = FALSE)
    }
    dup <- unique(as.character(tab[[id]][duplicated(tab[[id]])]))
    if (length(dup)) {
      stop("bindRich: x group '", group, "' ", slot, " lists identifier(s) ",
           .bind_quote(dup), " more than once.", call. = FALSE)
    }
  }
  absent <- "unavailable"
  if (has_tested) {
    testing <- enrichmentInfo(obj)$testing
    zoc <- if (is.list(testing)) testing$zero_overlap_count
    if (is.numeric(zoc) && length(zoc) == 1L && !is.na(zoc)) {
      absent <- if (zoc > 0) "untested_or_zero_overlap" else "not_tested"
    }
  }
  list(result = result, tested = tested, absent = absent)
}

#' Build one value column of the bindRich() table from the groups' own records
#' @return a vector (or list) of length n, typed as at first appearance
#' @keywords internal
#' @noRd
.bind_column <- function(col, records, groups, group_of, res_pos, tst_pos, n) {
  sources <- list()
  for (i in seq_along(records)) {
    r <- records[[i]]
    if (is.null(r)) next
    if (col %in% names(r$result)) {
      sources[[length(sources) + 1L]] <- list(group = i, from = "result",
                                              v = r$result[[col]])
    }
    if (!is.null(r$tested) && col %in% names(r$tested)) {
      sources[[length(sources) + 1L]] <- list(group = i, from = "tested",
                                              v = r$tested[[col]])
    }
  }
  types <- vapply(sources, function(s) .bind_type(s$v), character(1))
  numeric_family <- all(types %in% c("integer", "double"))
  if (!numeric_family && length(unique(types)) > 1L) {
    src_groups <- vapply(sources, `[[`, integer(1), "group")
    detail <- vapply(unique(src_groups), function(g) {
      paste0("'", groups[[g]], "' ",
             paste(unique(types[src_groups == g]), collapse = "/"))
    }, character(1))
    stop("bindRich: column '", col, "' has incompatible types across groups: ",
         paste(detail, collapse = ", "), ".", call. = FALSE)
  }

  first <- sources[[1L]]$v
  type <- types[[1L]]
  if (numeric_family && any(types == "double")) first <- as.double(first)
  out <- if (type == "list") {
    rep(list(NA_character_), n)
  } else if (type == "factor") {
    rep(NA_character_, n)
  } else {
    first[rep(NA_integer_, n)]
  }
  for (s in sources) {
    pos <- if (s$from == "result") res_pos else tst_pos
    take <- group_of == s$group & !is.na(pos)
    if (!any(take)) next
    v <- if (type == "factor") as.character(s$v) else s$v
    out[take] <- v[pos[take]]
  }
  if (type == "factor") {
    lv <- unique(unlist(lapply(sources, function(s) levels(s$v)),
                        use.names = FALSE))
    out <- factor(out, levels = lv)
  }
  out
}

#' Type signature of a column for the bindRich() type check
#' @keywords internal
#' @noRd
.bind_type <- function(v) {
  if (is.factor(v)) return("factor")
  if (is.list(v)) return("list")
  if (is.object(v)) return(paste(class(v), collapse = "/"))
  typeof(v)
}

#' @keywords internal
#' @noRd
.bind_quote <- function(x) paste0("'", x, "'", collapse = ", ")
