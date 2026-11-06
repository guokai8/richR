#' Convert gene identifiers with a recorded decision for every entry
#'
#' \code{mapGenes()} converts a gene vector or a ranking from one identifier
#' type to another and records what happened to every input entry. Unlike
#' \code{\link{idconvert}()}, which keeps the first candidate of a one-to-many
#' identifier without saying so, \code{mapGenes()} never selects a mapping
#' silently: each entry receives exactly one decision, and the decision table
#' travels with the result as its mapping record (read it with
#' \code{\link{mappingAudit}()}).
#'
#' @section Decisions:
#' Candidates of an entry are the distinct targets of the source pairs whose
#' \code{from} value equals the entry, in source order; a pair whose target is
#' missing or blank is not a candidate, and a repeated pair counts once. Each
#' entry leaves at the first step that decides it:
#' \enumerate{
#'   \item \code{invalid}: missing or blank identifier (for a ranking also a
#'     missing, \code{NaN} or infinite score);
#'   \item \code{unmapped}: no candidate;
#'   \item \code{ambiguous}: several candidates and \code{multiple = "drop"}
#'     (with \code{multiple = "first"} the first candidate in source order is
#'     kept and the entry continues; its \code{reason} names the candidate count);
#'   \item collision step over the entries still kept, grouped by target. For a
#'     gene vector the first entry of a target is \code{mapped} and later ones
#'     \code{merged}. For a ranking, entries sharing a target (including repeats
#'     of one source identifier) are kept once (\code{mapped}, later
#'     \code{merged}) when their scores are identical, and are all dropped with
#'     decision \code{conflict} when their scores differ. A score is never used
#'     to choose among entries.
#' }
#' The returned vector holds each kept target once, in the order of its first
#' kept entry; for a ranking its values are the input scores unchanged. Under
#' \code{multiple = "drop"} the result does not depend on the row order of the
#' source; under \code{"first"} the row order decides which candidate is kept.
#'
#' @section Mapping record:
#' The record is the attribute \code{"richR_mapping"} of the returned vector, a
#' plain list (named lists and atomic vectors only, no \code{class} or
#' \code{row.names}). It retains the conversion decisions and source facts
#' without imposing data-frame attributes:
#' \describe{
#'   \item{\code{decisions}}{the decision table as a named list of columns, one
#'     element per input entry in input order: \code{position}, \code{input_id},
#'     \code{candidates} (a list of character vectors, \code{character(0)} when
#'     none), \code{target} (\code{NA} when dropped), \code{decision},
#'     \code{reason}, and for a ranking \code{score} (the input score as given).}
#'   \item{\code{source}}{\code{kind} (\code{"table"} or \code{"OrgDb"}),
#'     \code{from}, \code{to}, \code{multiple}, \code{richR_version}; for an
#'     OrgDb source \code{package}, \code{version} and \code{metadata}
#'     (\code{list(name, value)} from \code{AnnotationDbi::metadata()}, rows in
#'     its order); for a table source \code{n_rows} and \code{fingerprint}.}
#' }
#' The table fingerprint is \code{list(algorithm = "sha256", canonicalization =
#' "richR-idmap-v1", hash, n_rows)}: the SHA-256 of the bytes of the text made
#' of the line \code{richR-idmap-v1} followed by one line \code{<from>\\t<to>}
#' per table row in row order (repeats kept), each line ended by LF, each field
#' UTF-8 with backslash, TAB, LF and CR written \code{\\\\}, \code{\\t},
#' \code{\\n}, \code{\\r} and a missing value written \code{\\N}.
#'
#' The record is lost by operations that drop attributes (\code{[},
#' \code{unique()}, \code{c()}, \code{sort()}). It is kept by \code{[<-};
#' \code{mappingAudit()} then reads the original conversion record, not a new
#' conversion of the modified vector. Existing enrichment functions do not
#' propagate this conversion record into their analysis metadata.
#'
#' @param x a character vector of identifiers, or a named numeric vector (a
#'   ranking, names are the identifiers).
#' @param mapping the mapping source: a data frame with character or factor
#'   columns named by \code{from} and \code{to}, an \code{OrgDb} object, or a
#'   species name listed by \code{\link{showData}()} (its OrgDb package must be
#'   installed).
#' @param from,to the source and target identifier types: column names of a
#'   table source, key types of an OrgDb source.
#' @param multiple how an entry with several candidates is resolved:
#'   \code{"drop"} (default, decision \code{ambiguous}) or \code{"first"} (the
#'   first candidate in source order).
#' @return a character vector of target identifiers, or a named numeric vector
#'   for a ranking, carrying the mapping record. When every entry is dropped
#'   the vector is empty and a warning states the count per decision.
#' @seealso \code{\link{mappingAudit}()}, \code{\link{idconvert}()}
#' @examples
#' tab <- data.frame(from = c("a", "b", "b", "c"), to = c("A", "B1", "B2", "A"))
#' v <- mapGenes(c("a", "b", "c", "x"), tab, "from", "to")
#' v
#' mappingAudit(v)
#' @export
mapGenes <- function(x, mapping, from, to, multiple = c("drop", "first")) {
  ranking <- .map_check_x(x)
  .map_check_type(from, "from")
  .map_check_type(to, "to")
  if (identical(from, to)) {
    stop(sprintf("mapGenes(): `from` and `to` are both \"%s\".", from), call. = FALSE)
  }
  if (identical(multiple, c("drop", "first"))) multiple <- "drop"
  if (!is.character(multiple) || length(multiple) != 1L || is.na(multiple) ||
      !multiple %in% c("drop", "first")) {
    stop(sprintf("mapGenes(): `multiple` must be \"drop\" or \"first\", not %s.",
                 paste(deparse(multiple), collapse = " ")), call. = FALSE)
  }
  ids <- if (ranking) names(x) else as.vector(x)
  ids <- as.character(ids)
  scores <- if (ranking) as.double(x) else NULL
  blank <- is.na(ids) | !nzchar(trimws(ids))
  invalid <- blank
  if (ranking) invalid <- invalid | !is.finite(scores)

  source <- .map_source(mapping, from, to, multiple)
  candidates <- source$candidates(unique(ids[!blank]))
  cand <- rep(list(character(0)), length(ids))
  found <- !blank & ids %in% names(candidates)
  cand[found] <- unname(candidates[ids[found]])
  n_cand <- lengths(cand)

  n <- length(ids)
  decision <- rep(NA_character_, n)
  reason <- rep(NA_character_, n)
  target <- rep(NA_character_, n)
  decision[invalid] <- "invalid"
  reason[blank] <- "missing or blank identifier"
  if (ranking) reason[!blank & invalid] <- "missing or non-finite score"
  todo <- is.na(decision)
  decision[todo & n_cand == 0L] <- "unmapped"
  reason[todo & n_cand == 0L] <- "no candidate in the mapping source"
  todo <- is.na(decision)
  amb <- todo & n_cand > 1L
  if (multiple == "drop") {
    decision[amb] <- "ambiguous"
    reason[amb] <- sprintf("%d candidates; multiple = \"drop\"", n_cand[amb])
  }
  kept <- which(is.na(decision))
  target[kept] <- vapply(cand[kept], `[[`, character(1), 1L)
  count_note <- ifelse(n_cand > 1L, sprintf("first of %d candidates", n_cand),
                       "1 candidate")
  # Collision step: entries still kept, grouped by target, in input order.
  groups <- split(kept, factor(target[kept], levels = unique(target[kept])))
  for (g in groups) {
    first <- g[[1L]]
    if (ranking && length(unique(scores[g])) > 1L) {
      decision[g] <- "conflict"
      reason[g] <- sprintf("%s; target %s shared by entries %s with different scores",
                           count_note[g], target[g], paste(g, collapse = ", "))
      target[g] <- NA_character_
      next
    }
    decision[first] <- "mapped"
    reason[first] <- count_note[first]
    if (length(g) > 1L) {
      rest <- g[-1L]
      decision[rest] <- "merged"
      reason[rest] <- sprintf("%s; same target as entry %d", count_note[rest], first)
    }
  }

  mapped <- which(decision == "mapped")
  out <- if (ranking) {
    stats::setNames(x[mapped], target[mapped])
  } else {
    target[mapped]
  }
  attributes(out) <- if (ranking) list(names = target[mapped]) else NULL
  decisions <- list(position = seq_len(n), input_id = ids, candidates = cand,
                    target = target, decision = decision, reason = reason)
  if (ranking) decisions$score <- scores
  attr(out, "richR_mapping") <- list(decisions = decisions, source = source$record)
  if (n > 0L && length(mapped) == 0L) {
    order <- c("invalid", "unmapped", "ambiguous", "conflict")
    counts <- table(factor(decision, levels = order))
    counts <- counts[counts > 0L]
    warning(sprintf("mapGenes(): every entry was dropped (%s).",
                    paste(names(counts), as.integer(counts), sep = ": ",
                          collapse = ", ")), call. = FALSE)
  }
  out
}

## TRUE for a ranking, FALSE for a gene vector; stops on anything else.
.map_check_x <- function(x) {
  what <- if (is.factor(x)) "factor" else if (is.list(x)) "list" else typeof(x)
  if (!is.null(dim(x)) || is.factor(x) || !(is.character(x) || is.numeric(x))) {
    stop(sprintf(paste("mapGenes(): `x` must be a character vector or a named",
                       "numeric vector, not a %s."), what), call. = FALSE)
  }
  if (length(x) == 0L) {
    stop(sprintf("mapGenes(): `x` is empty (a %s vector of length 0).", typeof(x)),
         call. = FALSE)
  }
  if (is.numeric(x) && is.null(names(x))) {
    stop("mapGenes(): `x` is a numeric vector without names; a ranking needs gene names.",
         call. = FALSE)
  }
  is.numeric(x)
}

.map_check_type <- function(value, arg) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(trimws(value))) {
    stop(sprintf("mapGenes(): `%s` must be one identifier type (a non-blank string), not %s.",
                 arg, paste(deparse(value), collapse = " ")), call. = FALSE)
  }
}

## The mapping source: its record (identity) and a candidate lookup that
## returns, for unique valid ids, a named list of distinct candidate targets in
## source order.
.map_source <- function(mapping, from, to, multiple) {
  base <- list(from = from, to = to, multiple = multiple,
               richR_version = .richr_version())
  if (is.data.frame(mapping)) return(.map_table_source(mapping, from, to, base))
  if (is.character(mapping) && length(mapping) == 1L && !is.na(mapping)) {
    package <- .getdb(mapping)
    if (is.null(package)) {
      stop(sprintf("mapGenes(): `mapping` = \"%s\" is not a species listed by showData().",
                   mapping), call. = FALSE)
    }
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(sprintf(paste("mapGenes(): `mapping` = \"%s\" needs the package \"%s\",",
                         "which is not installed."), mapping, package), call. = FALSE)
    }
    mapping <- getExportedValue(package, package)
  }
  if (methods::is(mapping, "OrgDb")) return(.map_orgdb_source(mapping, from, to, base))
  stop(sprintf(paste("mapGenes(): `mapping` must be a data frame, an OrgDb object or a",
                     "species name listed by showData(), not a %s."),
               class(mapping)[[1L]]), call. = FALSE)
}

.map_table_source <- function(mapping, from, to, base) {
  for (arg in c("from", "to")) {
    col <- if (arg == "from") from else to
    if (!col %in% names(mapping)) {
      stop(sprintf("mapGenes(): `%s` = \"%s\" is not a column of `mapping` (columns: %s).",
                   arg, col, paste(names(mapping), collapse = ", ")), call. = FALSE)
    }
    v <- mapping[[col]]
    if (!(is.character(v) || is.factor(v))) {
      stop(sprintf(paste("mapGenes(): column `%s` (`%s` = \"%s\") of `mapping` is %s;",
                         "identifier columns must be character or factor."),
                   col, arg, col, class(v)[[1L]]), call. = FALSE)
    }
  }
  f <- as.character(mapping[[from]])
  t <- as.character(mapping[[to]])
  ok <- !is.na(f) & !is.na(t) & nzchar(trimws(t))
  pairs <- split(t[ok], factor(f[ok], levels = unique(f[ok])))
  lookup <- lapply(pairs, function(v) unique(unname(v)))
  record <- c(list(kind = "table"), base,
              list(n_rows = nrow(mapping),
                   fingerprint = list(algorithm = "sha256",
                                      canonicalization = "richR-idmap-v1",
                                      hash = .idmap_hash(f, t),
                                      n_rows = nrow(mapping))))
  list(record = record, candidates = function(ids) lookup[intersect(ids, names(lookup))])
}

## SHA-256 of the richR-idmap-v1 canonical text (see ?mapGenes).
.idmap_hash <- function(from, to) {
  esc <- function(v) {
    v <- enc2utf8(v)
    v <- gsub("\\", "\\\\", v, fixed = TRUE)
    v <- gsub("\t", "\\t", v, fixed = TRUE)
    v <- gsub("\n", "\\n", v, fixed = TRUE)
    v <- gsub("\r", "\\r", v, fixed = TRUE)
    v[is.na(v)] <- "\\N"
    v
  }
  text <- paste0("richR-idmap-v1\n", paste0(esc(from), "\t", esc(to), "\n", collapse = ""))
  digest(text, algo = "sha256", serialize = FALSE)
}

.map_orgdb_source <- function(db, from, to, base) {
  if (!from %in% AnnotationDbi::keytypes(db)) {
    stop(sprintf("mapGenes(): `from` = \"%s\" is not a key type of the OrgDb source.", from),
         call. = FALSE)
  }
  if (!to %in% AnnotationDbi::columns(db)) {
    stop(sprintf("mapGenes(): `to` = \"%s\" is not a key type of the OrgDb source.", to),
         call. = FALSE)
  }
  package <- db$packageName
  if (!is.character(package) || length(package) != 1L) package <- NA_character_
  md <- AnnotationDbi::metadata(db)
  record <- c(list(kind = "OrgDb"), base,
              list(package = package,
                   version = .package_version_or_na(package),
                   metadata = list(name = as.character(md$name),
                                   value = as.character(md$value))))
  candidates <- function(ids) {
    ids <- ids[ids %in% AnnotationDbi::keys(db, keytype = from)]
    if (!length(ids)) return(list())
    m <- suppressMessages(AnnotationDbi::mapIds(db, keys = ids, column = to,
                                                keytype = from, multiVals = "list"))
    lapply(m, function(v) {
      v <- unname(as.character(v))
      unique(v[!is.na(v) & nzchar(trimws(v))])
    })
  }
  list(record = record, candidates = candidates)
}
