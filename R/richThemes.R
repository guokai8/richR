#' Group reported terms into themes by the overlap of their evidence genes
#'
#' Groups the terms that one or several enrichment results report into
#' themes of related terms, by the kappa similarity of the genes that support
#' each term, and keeps everything needed to go back from a theme to its
#' terms, their rows, their genes and the input statistics. A theme is a
#' summary of the listed terms, not a statistical test.
#'
#' @section Interpretation:
#' \itemize{
#'   \item Theme grouping summarizes the listed terms and is \strong{not a
#'     statistical test}. No theme-level p-value exists: no element of the
#'     object carries a p-value, adjusted p-value, NES, enrichment score or
#'     any value computed from them, and the representative's p-value is not
#'     a p-value of the theme.
#'   \item Related or overlapping gene sets (for example the same pathway in
#'     two databases, or a parent and a child GO term) are \strong{not
#'     independent replications}: they share genes by construction, so the
#'     number of terms in a theme is not a measure of evidence.
#'   \item Evidence genes are the ORA overlap genes or the GSEA leading-edge
#'     genes of a term. They are the genes behind the term's result,
#'     \strong{not drivers} shown to cause anything.
#' }
#'
#' @section Grouping rule:
#' Every input follows the group rules of \code{\link{bindRich}}: a single
#' result is treated as \code{list(x)} and named \code{Group 1}; an unnamed
#' list is named \code{Group 1}, \code{Group 2}, ...; \code{NULL} entries are
#' failed groups; ORA and GSEA results cannot be mixed. Identifiers are
#' \code{Annot} (ORA) or \code{pathway} (GSEA); \code{Term} is display text
#' only.
#'
#' \describe{
#'   \item{Grouping set}{the identifiers whose \code{bindRich()} state is
#'     \code{"reported"} in at least one group. Terms that are only tested,
#'     not tested, \code{"untested_or_zero_overlap"} or \code{"unavailable"}
#'     are never grouped.}
#'   \item{Evidence genes}{of a term in a group, read from that group's own
#'     record: its \code{result(x)} row when the term is reported there, its
#'     tested-table row when it is tested but not reported there, nothing
#'     otherwise. ORA: the \code{GeneID} cell split on the object's
#'     \code{sep} (\code{","} when the slot is empty). GSEA: the
#'     \code{leadingEdge} cell (a list element, or text split on \code{sep}).
#'     Missing and empty strings are not genes; duplicates count once.}
#'   \item{Similarity genes}{of a grouped term: the union of its evidence
#'     genes over the groups that report it. They pool the supplied groups,
#'     so \strong{themes can change when a group is added or removed}.}
#'   \item{Background}{the union of the similarity genes of every grouped
#'     term; its size is \code{settings$background_size}.}
#'   \item{Similarity}{Cohen's kappa of two terms' similarity genes against
#'     the background (the package's kappa statistic, as used by
#'     \code{\link{richCluster}}). No p-value, adjusted p-value, NES or
#'     enrichment score is used anywhere in grouping.}
#'   \item{Clustering}{average-linkage hierarchical clustering of
#'     \code{1 - kappa}, cut at height \code{1 - cutoff}. The cut is
#'     inclusive: two clusters join when their average kappa is \strong{at
#'     least} \code{cutoff}, whereas \code{richCluster()} keeps pairs with
#'     kappa strictly greater than its cutoff, so the same number does not
#'     mean the same threshold in the two functions. Every grouped term
#'     belongs to exactly one theme; a term joined to no other is a singleton
#'     theme; no theme is removed by size or by any p-value.}
#'   \item{Degenerate kappa}{kappa is \code{0} when two terms share no gene,
#'     and also when the chance agreement is 1 — for example when every
#'     grouped term has the same genes and the background holds no other gene
#'     (kappa cannot be computed there). Such terms stay singletons; the
#'     small \code{background_size} shows why. A term without evidence genes
#'     is a singleton with no \code{genes} row.}
#'   \item{Representative}{the member with the highest mean kappa to the
#'     other members; means within \code{1e-12} of the highest are ties.
#'     Ties and singletons resolve to the first identifier in radix (C-locale)
#'     order.}
#'   \item{Label}{the representative's \code{Term} in the radix-first group
#'     name among the groups that report it with a non-missing, non-blank
#'     \code{Term} (\code{representative_term}); when there is none, the
#'     identifier. Labels are only display text: edit one with
#'     \code{th$themes$label[th$themes$theme_id == "theme_001"] <- "my label"};
#'     every link between tables uses \code{theme_id}.}
#' }
#'
#' Themes, theme ids, representatives, labels and the \code{similarity}
#' table do not depend on the order of the groups or of the rows within a
#' group. The \code{group} factor levels, the order of rows by group, the
#' order of the \code{members} columns after \code{state} and
#' \code{execution} and the levels of factor columns inherited from
#' \code{bindRich()}, and \code{settings$groups} follow the input order.
#'
#' @section Value:
#' A list of class \code{"richThemes"} with six elements:
#' \describe{
#'   \item{\code{themes}}{one row per theme: \code{theme_id} (\code{"theme_"}
#'     and the theme's position, zero-padded to at least three digits),
#'     \code{label}, \code{representative} (identifier),
#'     \code{representative_term} (may be \code{NA}), \code{n_terms},
#'     \code{n_groups_reported} (groups where at least one member is
#'     reported), \code{n_similarity_genes} (size of the union of the
#'     members' similarity genes, so genes of reported records only;
#'     \code{\link{themeGenes}} may list more, because it also shows
#'     tested-not-reported records) and \code{mean_within_kappa} (mean
#'     pairwise kappa among members, \code{NA} for a singleton). Themes are
#'     ordered by \code{n_terms} descending, then by representative.}
#'   \item{\code{members}}{one row per member and group: \code{theme_id},
#'     \code{is_representative}, then every \code{bindRich(x)} column in its
#'     order, with its values and classes unchanged. Rows are ordered by
#'     theme, identifier and group, and are exactly the \code{bindRich(x)}
#'     rows of the members.}
#'   \item{\code{genes}}{one row per member, group and evidence gene, for the
#'     groups where the member is reported or tested but not reported:
#'     \code{theme_id}, \code{group}, the identifier, \code{gene},
#'     \code{evidence_type} (\code{"ora_overlap"} or
#'     \code{"gsea_leading_edge"}) and \code{state}.}
#'   \item{\code{similarity}}{one row per unordered pair of members of one
#'     theme: \code{theme_id}, \code{identifier_a}, \code{identifier_b}
#'     (radix order), \code{kappa} and \code{n_shared} (similarity genes the
#'     two share).}
#'   \item{\code{ranking}}{GSEA: one row per group and ranked gene, the
#'     group's ranked input in stored order: \code{group}, \code{gene},
#'     \code{rank_score}. Zero rows for ORA.}
#'   \item{\code{settings}}{\code{method} (\code{"kappa"}), \code{distance}
#'     (\code{"1 - kappa"}), \code{linkage} (\code{"average"}),
#'     \code{cutoff}, \code{grouped_on}, \code{similarity_genes},
#'     \code{background} (the rules above, as text), \code{background_size},
#'     \code{representative}, \code{label}, \code{groups} (input order),
#'     \code{sources} (each group's \code{\link{enrichmentInfo}} record,
#'     \code{"failed"} for a \code{NULL} group) and \code{richR_version}.}
#' }
#' When nothing is reported in any group, the five tables have zero rows with
#' their columns and classes, and \code{settings} is filled. Every column
#' \code{richThemes()} creates is character, integer, double, logical or
#' factor; \code{members} keeps the classes of \code{bindRich()}.
#'
#' @param x a \code{richResult}, a \code{GSEAResult}, or a list of them (with
#'   \code{NULL} for failed groups), as accepted by \code{\link{bindRich}}.
#' @param cutoff the kappa threshold, one finite number in (0, 1]: clusters
#'   join when their average kappa is at least \code{cutoff}. Default 0.5.
#' @return a \code{richThemes} object, described in the Value section.
#' @seealso \code{\link{themeGenes}}, \code{\link{bindRich}},
#'   \code{\link{richCluster}}
#' @examples
#' \dontrun{
#' res <- batchEnrich(list(S1 = gene1, S2 = gene2), annot = hsako)
#' th <- richThemes(res$results)
#' th$themes
#' th$themes$label[1] <- "my label"
#' themeGenes(th, "theme_001")
#' }
#' @export
richThemes <- function(x, cutoff = 0.5) {
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) ||
      cutoff <= 0 || cutoff > 1) {
    stop("richThemes: cutoff must be one finite number in (0, 1].",
         call. = FALSE)
  }
  if (methods::is(x, "richResult") || methods::is(x, "GSEAResult")) {
    x <- list(x)
  }
  bound <- bindRich(x)
  if (is.null(names(x))) names(x) <- paste("Group", seq_along(x))
  groups <- names(x)
  gsea <- any(vapply(x, methods::is, logical(1), "GSEAResult"))
  id <- if (gsea) "pathway" else "Annot"
  added <- c("theme_id", "is_representative")
  clash <- intersect(names(bound), added)
  if (length(clash)) {
    stop("richThemes: result column(s) ", .bind_quote(clash),
         " clash with the columns richThemes() adds.", call. = FALSE)
  }

  ids <- as.character(bound[[id]])
  state <- bound$state
  ev <- .theme_evidence(bound, x, groups, gsea)

  # Grouping set, similarity genes and background (radix order throughout,
  # so nothing depends on the input order or the session locale).
  grouped <- sort(unique(ids[state %in% "reported"]), method = "radix")
  rep_rows <- which(state %in% "reported")
  sim_genes <- lapply(split(rep_rows, factor(ids[rep_rows], levels = grouped)),
                      function(r) {
                        g <- as.character(unlist(ev[r], use.names = FALSE))
                        sort(unique(g), method = "radix")
                      })
  names(sim_genes) <- grouped
  background <- sort(unique(as.character(unlist(sim_genes, use.names = FALSE))),
                     method = "radix")
  kap <- .theme_kappa(sim_genes, background)

  n <- length(grouped)
  cluster <- if (n >= 2L) {
    stats::cutree(stats::hclust(stats::as.dist(1 - kap), method = "average"),
                  h = 1 - cutoff)
  } else {
    rep(1L, n)
  }
  mem <- unname(split(grouped, cluster))
  rep_of <- vapply(mem, .theme_representative, character(1), kap = kap)
  ord <- order(-lengths(mem), rep_of, method = "radix")
  mem <- mem[ord]
  rep_of <- rep_of[ord]
  n_themes <- length(mem)
  width <- max(3L, nchar(as.character(n_themes)))
  theme_ids <- sprintf(paste0("theme_%0", width, "d"), seq_len(n_themes))

  rep_term <- vapply(rep_of, .theme_rep_term, character(1), ids = ids,
                     state = state, bound = bound)
  rep_term <- unname(rep_term)
  label <- ifelse(is.na(rep_term), rep_of, rep_term)
  themes <- .theme_frame(list(
    theme_id = theme_ids,
    label = as.character(label),
    representative = as.character(rep_of),
    representative_term = as.character(rep_term),
    n_terms = as.integer(lengths(mem)),
    n_groups_reported = vapply(mem, function(m) {
      length(unique(as.character(bound$group[ids %in% m & state %in% "reported"])))
    }, integer(1)),
    n_similarity_genes = vapply(mem, function(m) {
      length(unique(unlist(sim_genes[m], use.names = FALSE)))
    }, integer(1)),
    mean_within_kappa = vapply(mem, function(m) {
      if (length(m) < 2L) return(NA_real_)
      k <- kap[m, m, drop = FALSE]
      mean(k[upper.tri(k)])
    }, numeric(1))
  ))

  # members: the bindRich() rows of every member in every group.
  member_ids <- unlist(mem, use.names = FALSE)
  theme_index <- rep(seq_len(n_themes), lengths(mem))
  theme_of <- function(v) theme_index[match(v, member_ids)]
  sel <- which(ids %in% grouped)
  sel <- sel[order(theme_of(ids[sel]), match(ids[sel], grouped),
                   as.integer(bound$group[sel]))]
  members <- .theme_frame(c(
    list(theme_id = theme_ids[theme_of(ids[sel])],
         is_representative = ids[sel] %in% rep_of),
    as.list(bound[sel, , drop = FALSE])
  ))

  # genes: each member's evidence genes in each group where it is reported
  # or tested but not reported, in the members' row order.
  gsel <- sel[state[sel] %in% c("reported", "tested_not_reported")]
  len <- lengths(ev[gsel])
  genes <- .theme_frame(stats::setNames(list(
    theme_ids[rep(theme_of(ids[gsel]), len)],
    rep(bound$group[gsel], len),
    rep(ids[gsel], len),
    as.character(unlist(ev[gsel], use.names = FALSE)),
    rep(if (gsea) "gsea_leading_edge" else "ora_overlap", sum(len)),
    rep(state[gsel], len)
  ), c("theme_id", "group", id, "gene", "evidence_type", "state")))

  similarity <- .theme_similarity(mem, theme_ids, kap, sim_genes)
  ranking <- .theme_ranking(x, groups, gsea)

  settings <- list(
    method = "kappa",
    distance = "1 - kappa",
    linkage = "average",
    cutoff = cutoff,
    grouped_on = paste("identifiers whose bindRich() state is 'reported'",
                       "in at least one group"),
    similarity_genes = paste(
      "a term's evidence genes unioned over the groups that report it;",
      if (gsea) "evidence genes are the leadingEdge genes (GSEA)" else
        "evidence genes are the GeneID overlap genes (ORA)",
      "of the group's own record"),
    background = "the union of the similarity genes of every grouped term",
    background_size = length(background),
    representative = paste(
      "the member with the highest mean kappa to the other members; means",
      "within 1e-12 of the highest are ties; ties and singletons resolve to",
      "the first identifier in radix order"),
    label = paste(
      "the representative's Term in the radix-first group that reports it",
      "with a non-missing, non-blank Term; otherwise the identifier"),
    groups = groups,
    sources = lapply(x, function(o) if (is.null(o)) "failed" else enrichmentInfo(o)),
    richR_version = as.character(utils::packageVersion("richR"))
  )

  structure(list(themes = themes, members = members, genes = genes,
                 similarity = similarity, ranking = ranking,
                 settings = settings),
            class = "richThemes")
}

#' A data.frame from a named list of equal-length columns, with compact row
#' names and list columns kept as they are
#' @keywords internal
#' @noRd
.theme_frame <- function(cols) {
  n <- if (length(cols)) length(cols[[1L]]) else 0L
  structure(cols, names = names(cols), row.names = .set_row_names(n),
            class = "data.frame")
}

#' The separator a result object records, or "," when it records none
#' @keywords internal
#' @noRd
.theme_sep <- function(obj) {
  s <- obj@sep
  if (length(s) == 1L && !is.na(s) && nzchar(s)) s else ","
}

#' Evidence genes of every bindRich() row that is reported or tested but not
#' reported, read from the row's own record (bindRich() takes each value from
#' the group's result or tested table); character(0) for every other row
#' @keywords internal
#' @noRd
.theme_evidence <- function(bound, x, groups, gsea) {
  out <- rep(list(character(0)), nrow(bound))
  col <- if (gsea) "leadingEdge" else "GeneID"
  if (!col %in% names(bound)) return(out)
  v <- bound[[col]]
  seps <- vapply(x, function(o) if (is.null(o)) "," else .theme_sep(o), character(1))
  gi <- as.integer(bound$group)
  for (r in which(bound$state %in% c("reported", "tested_not_reported"))) {
    if (is.list(v)) {
      g <- as.character(v[[r]])
    } else {
      cell <- as.character(v[[r]])
      g <- if (is.na(cell)) character(0) else
        strsplit(cell, seps[[gi[[r]]]], fixed = TRUE)[[1L]]
    }
    g <- g[!is.na(g) & nzchar(g)]
    out[[r]] <- unique(g)
  }
  out
}

#' Cohen's kappa between every pair of gene sets against a background,
#' computed from the 2 x 2 table exactly as .kappa() does, on gene vectors
#' @keywords internal
#' @noRd
.theme_kappa <- function(sets, background) {
  n <- length(sets)
  N <- length(background)
  inc <- matrix(0, N, n)
  for (i in seq_len(n)) inc[match(sets[[i]], background), i] <- 1
  n11 <- crossprod(inc)
  n10 <- matrix(colSums(inc), n, n) - n11
  n01 <- t(n10)
  n00 <- N - n11 - n10 - n01
  oab <- (n11 + n00) / N
  aab <- ((n11 + n10) * (n11 + n01) + (n01 + n00) * (n10 + n00)) / (N * N)
  kap <- (oab - aab) / (1 - aab)
  kap[n11 == 0 | aab == 1] <- 0
  diag(kap) <- 1
  dimnames(kap) <- list(names(sets), names(sets))
  kap
}

#' The member with the highest mean kappa to the others (ties within 1e-12
#' resolve to the first identifier; members arrive in radix order)
#' @keywords internal
#' @noRd
.theme_representative <- function(m, kap) {
  if (length(m) == 1L) return(m)
  k <- kap[m, m, drop = FALSE]
  diag(k) <- NA
  mk <- rowMeans(k, na.rm = TRUE)
  m[mk >= max(mk) - 1e-12][[1L]]
}

#' The representative's Term in the radix-first group that reports it with a
#' non-missing, non-blank Term, or NA
#' @keywords internal
#' @noRd
.theme_rep_term <- function(r, ids, state, bound) {
  rows <- which(ids == r & state %in% "reported")
  term <- as.character(bound$Term[rows])
  grp <- as.character(bound$group[rows])
  ok <- !is.na(term) & nzchar(trimws(term))
  if (!any(ok)) return(NA_character_)
  term[ok][order(grp[ok], method = "radix")][[1L]]
}

#' One row per unordered pair of members of one theme
#' @keywords internal
#' @noRd
.theme_similarity <- function(mem, theme_ids, kap, sim_genes) {
  parts <- lapply(seq_along(mem), function(t) {
    m <- mem[[t]]
    if (length(m) < 2L) return(NULL)
    p <- which(upper.tri(diag(length(m))), arr.ind = TRUE)
    p <- p[order(p[, 1L], p[, 2L]), , drop = FALSE]
    a <- m[p[, 1L]]
    b <- m[p[, 2L]]
    list(theme_id = rep(theme_ids[[t]], length(a)), identifier_a = a,
         identifier_b = b, kappa = kap[cbind(a, b)],
         n_shared = vapply(seq_along(a), function(j) {
           length(intersect(sim_genes[[a[[j]]]], sim_genes[[b[[j]]]]))
         }, integer(1)))
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  .theme_frame(list(
    theme_id = as.character(unlist(lapply(parts, `[[`, "theme_id"))),
    identifier_a = as.character(unlist(lapply(parts, `[[`, "identifier_a"))),
    identifier_b = as.character(unlist(lapply(parts, `[[`, "identifier_b"))),
    kappa = as.numeric(unlist(lapply(parts, `[[`, "kappa"))),
    n_shared = as.integer(unlist(lapply(parts, `[[`, "n_shared")))
  ))
}

#' Every GSEA group's ranked input in stored order; zero rows for ORA
#' @keywords internal
#' @noRd
.theme_ranking <- function(x, groups, gsea) {
  keep <- if (gsea) which(!vapply(x, is.null, logical(1))) else integer(0)
  inputs <- lapply(x[keep], function(o) o@input)
  gene <- lapply(inputs, function(v) {
    if (is.null(names(v))) rep(NA_character_, length(v)) else names(v)
  })
  .theme_frame(list(
    group = factor(rep(groups[keep], lengths(inputs)), levels = groups),
    gene = as.character(unlist(gene, use.names = FALSE)),
    rank_score = as.numeric(unlist(lapply(inputs, unname), use.names = FALSE))
  ))
}

#' Show the evidence genes of themes per group, with input statistics
#'
#' Drills down from themes made by \code{\link{richThemes}} to their genes:
#' one row per theme, evidence gene and group, saying which members list the
#' gene in that group, what the theme's status is in that group, and the
#' gene's statistics when they are supplied. Statistics are only joined; no
#' similarity, overlap test, score or difference between groups is computed.
#'
#' The genes of a theme are every evidence gene of any member in any group:
#' ORA overlap or GSEA leading-edge genes of reported and of tested but not
#' reported records (so there may be more genes than the theme's
#' \code{n_similarity_genes}, which counts reported records only). Evidence
#' genes are the genes behind a term's result, \strong{not drivers} shown to
#' cause anything. Grouping terms into themes is \strong{not a statistical
#' test} and \strong{no theme-level p-value exists}; related or overlapping
#' gene sets are \strong{not independent replications}, so the number of
#' members listing a gene is not a measure of evidence. See
#' \code{\link{richThemes}}.
#'
#' @section Value:
#' A data.frame with one row per theme, gene and group, ordered by theme (as
#' requested), gene (radix order) and group (input order), with the columns:
#' \describe{
#'   \item{\code{theme_id}, \code{label}}{the theme and its label as
#'     currently stored in \code{themes$themes}, so an edited label shows;}
#'   \item{\code{gene}}{the evidence gene;}
#'   \item{\code{group}}{a factor with the groups in input order;}
#'   \item{\code{evidence_members}}{the members whose record in this group
#'     lists the gene, comma-joined in radix order, \code{""} when none; and
#'     \code{n_evidence_members}, their number. A theme seen in two groups is
#'     supported by the same genes only where both groups list them;}
#'   \item{\code{group_status}}{the theme's status in the group, first match
#'     wins: \code{"failed"} (a \code{NULL} group), \code{"reported"} (a
#'     member is reported), \code{"tested_only"} (a member is tested but not
#'     reported), \code{"unavailable"} (a member's state cannot be known),
#'     \code{"undetermined"} (a member is untested or has zero overlap and
#'     the tested table cannot distinguish these), \code{"absent"} otherwise.
#'     Here \code{"unavailable"} remains separate; the \code{undetermined}
#'     count in \code{\link{richSensitivity}} also includes unavailable states;}
#'   \item{the statistic columns}{the further columns of \code{stats}, joined
#'     by exact gene string; for GSEA without \code{stats}, \code{rank_score},
#'     the group's ranked input; for ORA without \code{stats}, none;}
#'   \item{\code{stat_status}}{why a statistic is or is not there, first match
#'     wins: \code{"group_failed"} (a \code{NULL} group, even when
#'     \code{stats} has an entry for it), \code{"no_statistics_for_group"},
#'     \code{"joined"}, \code{"gene_not_in_statistics"}. A statistic that is
#'     not joined is \code{NA}, never \code{0}.}
#' }
#'
#' @param themes a \code{richThemes} object.
#' @param theme \code{NULL} for every theme, or one or more \code{theme_id}
#'   values to show; an empty vector is an error.
#' @param stats \code{NULL}, or a named list of data frames, one per group
#'   (names among the object's groups; groups may be left out). Each needs a
#'   character column \code{gene} listing each gene once and at least one
#'   further column; every entry has the same further columns, each of the
#'   same class in every entry (and, for a factor, with the same levels, so
#'   no value is coerced or lost when the groups are joined), and none may
#'   be named like an output column (\code{theme_id}, \code{label},
#'   \code{group}, \code{evidence_members}, \code{n_evidence_members},
#'   \code{group_status}, \code{stat_status}). \code{stats} that break these
#'   rules are refused with an error naming the offending entry or column.
#'   Gene identifiers are not mapped: use the identifiers of the results.
#' @return a data.frame, described in the Value section.
#' @seealso \code{\link{richThemes}}
#' @examples
#' \dontrun{
#' th <- richThemes(list(A = resA, B = resB))
#' fc <- list(A = data.frame(gene = names(fcA), log2FC = fcA),
#'            B = data.frame(gene = names(fcB), log2FC = fcB))
#' themeGenes(th, "theme_001", stats = fc)
#' }
#' @export
themeGenes <- function(themes, theme = NULL, stats = NULL) {
  if (!inherits(themes, "richThemes")) {
    stop("themeGenes: themes must be a richThemes object, as returned by ",
         "richThemes().", call. = FALSE)
  }
  all_ids <- themes$themes$theme_id
  if (is.null(theme)) {
    theme <- all_ids
  } else {
    if (!length(theme)) {
      stop("themeGenes: theme must name at least one theme_id of themes, or be ",
           "NULL for every theme.", call. = FALSE)
    }
    bad <- setdiff(as.character(theme), all_ids)
    if (!is.character(theme) || length(bad)) {
      stop("themeGenes: theme ",
           .bind_quote(if (length(bad)) bad else as.character(theme)),
           " is not a theme_id of themes.", call. = FALSE)
    }
    theme <- unique(theme)
  }
  groups <- themes$settings$groups
  id <- names(themes$genes)[[3L]]
  scols <- character(0)
  template <- NULL
  if (!is.null(stats)) {
    scols <- .theme_check_stats(stats, groups)
    template <- stats[[1L]][0L, scols, drop = FALSE]
  } else if (identical(id, "pathway")) {
    scols <- "rank_score"
    rk <- themes$ranking
    stats <- lapply(stats::setNames(nm = unique(as.character(rk$group))),
                    function(g) {
                      sel <- as.character(rk$group) == g
                      data.frame(gene = rk$gene[sel],
                                 rank_score = rk$rank_score[sel],
                                 stringsAsFactors = FALSE)
                    })
    template <- data.frame(rank_score = numeric(0))
  }

  parts <- lapply(theme, .theme_gene_rows, themes = themes, groups = groups,
                  id = id, stats = stats, scols = scols, template = template)
  col <- function(name) do.call(c, lapply(parts, `[[`, name))
  out <- list(
    theme_id = as.character(col("theme_id")),
    label = as.character(col("label")),
    gene = as.character(col("gene")),
    group = factor(groups[as.integer(col("group"))], levels = groups),
    evidence_members = as.character(col("evidence_members")),
    n_evidence_members = as.integer(col("n_evidence_members")),
    group_status = as.character(col("group_status"))
  )
  for (cc in scols) {
    v <- col(cc)
    out[[cc]] <- if (is.null(v)) template[[cc]] else v
  }
  out$stat_status <- as.character(col("stat_status"))
  .theme_frame(out)
}

#' Validate themeGenes() statistics; returns the further column names
#' @keywords internal
#' @noRd
.theme_check_stats <- function(stats, groups) {
  nm <- names(stats)
  if (!is.list(stats) || is.data.frame(stats) || !length(stats) || is.null(nm) ||
      anyNA(nm) || any(!nzchar(trimws(nm)))) {
    stop("themeGenes: stats must be NULL or a named list of data frames, ",
         "one per group.", call. = FALSE)
  }
  dup <- unique(nm[duplicated(nm)])
  if (length(dup)) {
    stop("themeGenes: stats names group(s) ", .bind_quote(dup),
         " more than once.", call. = FALSE)
  }
  unknown <- setdiff(nm, groups)
  if (length(unknown)) {
    stop("themeGenes: stats names group(s) ", .bind_quote(unknown),
         " that themes does not hold (groups: ", .bind_quote(groups), ").",
         call. = FALSE)
  }
  reserved <- c("theme_id", "label", "group", "evidence_members",
                "n_evidence_members", "group_status", "stat_status")
  further <- list()
  for (g in nm) {
    d <- stats[[g]]
    if (!is.data.frame(d)) {
      stop("themeGenes: stats entry '", g, "' must be a data frame.",
           call. = FALSE)
    }
    if (!"gene" %in% names(d) || !is.character(d$gene)) {
      stop("themeGenes: stats entry '", g, "' needs a character column ",
           "'gene'.", call. = FALSE)
    }
    dg <- unique(d$gene[duplicated(d$gene)])
    if (length(dg)) {
      stop("themeGenes: stats entry '", g, "' lists gene(s) ", .bind_quote(dg),
           " more than once.", call. = FALSE)
    }
    f <- setdiff(names(d), "gene")
    if (!length(f)) {
      stop("themeGenes: stats entry '", g, "' has no column besides 'gene'.",
           call. = FALSE)
    }
    clash <- intersect(f, reserved)
    if (length(clash)) {
      stop("themeGenes: stats entry '", g, "' column(s) ", .bind_quote(clash),
           " clash with the columns themeGenes() returns.", call. = FALSE)
    }
    further[[g]] <- f
  }
  differ <- setdiff(unique(unlist(further, use.names = FALSE)),
                    Reduce(intersect, further))
  if (length(differ)) {
    stop("themeGenes: stats entries differ in their column(s) ",
         .bind_quote(differ), "; every entry needs the same columns.",
         call. = FALSE)
  }
  # Each statistic column is joined into one vector, so a class or factor
  # levels that differ between entries would coerce values or turn them NA.
  same <- "; every entry needs the same class and, for a factor, the same levels."
  for (cc in further[[1L]]) {
    cls <- vapply(stats, function(d) paste(class(d[[cc]]), collapse = "/"),
                  character(1))
    if (length(unique(cls)) > 1L) {
      stop("themeGenes: stats column ", .bind_quote(cc), " has a different class ",
           "in different entries (", paste0("'", nm, "': ", cls, collapse = "; "),
           ")", same, call. = FALSE)
    }
    lev <- lapply(stats, function(d) levels(d[[cc]]))
    off <- !vapply(lev, identical, logical(1), lev[[1L]])
    if (any(off)) {
      stop("themeGenes: stats column ", .bind_quote(cc), " has different ",
           "factor levels in entries ", .bind_quote(c(nm[[1L]], nm[off])), same,
           call. = FALSE)
    }
  }
  further[[1L]]
}

#' The themeGenes() rows of one theme, as a list of columns ordered by gene
#' (radix) and group (input order)
#' @keywords internal
#' @noRd
.theme_gene_rows <- function(t, themes, groups, id, stats, scols, template) {
  gm <- themes$genes[themes$genes$theme_id == t, , drop = FALSE]
  genes <- sort(unique(as.character(gm$gene)), method = "radix")
  ng <- length(groups)
  nu <- length(genes)
  mm <- themes$members[themes$members$theme_id == t, , drop = FALSE]
  mgroup <- as.character(mm$group)
  status <- vapply(groups, function(g) {
    st <- mm$state[mgroup == g]
    if (any(mm$execution[mgroup == g] %in% "failed")) return("failed")
    if (any(st %in% "reported")) return("reported")
    if (any(st %in% "tested_not_reported")) return("tested_only")
    if (any(st %in% "unavailable")) return("unavailable")
    if (any(st %in% "untested_or_zero_overlap")) return("undetermined")
    "absent"
  }, character(1), USE.NAMES = FALSE)

  # Built group by group, then reordered to gene, group.
  ev <- lapply(groups, function(g) {
    sel <- as.character(gm$group) == g
    split(as.character(gm[[id]][sel]), factor(gm$gene[sel], levels = genes))
  })
  em <- lapply(ev, vapply, function(v) {
    paste(sort(unique(v), method = "radix"), collapse = ",")
  }, character(1), USE.NAMES = FALSE)
  nem <- lapply(ev, vapply, function(v) length(unique(v)), integer(1),
                USE.NAMES = FALSE)
  sstat <- vector("list", ng)
  svals <- lapply(stats::setNames(nm = scols), function(cc) vector("list", ng))
  for (k in seq_len(ng)) {
    d <- if (!is.null(stats)) stats[[groups[[k]]]] else NULL
    failed <- status[[k]] == "failed"
    pos <- if (is.null(d) || failed) rep(NA_integer_, nu) else match(genes, d$gene)
    sstat[[k]] <- if (failed) {
      rep("group_failed", nu)
    } else if (is.null(d)) {
      rep("no_statistics_for_group", nu)
    } else {
      ifelse(is.na(pos), "gene_not_in_statistics", "joined")
    }
    for (cc in scols) {
      v <- template[[cc]][rep(NA_integer_, nu)]
      take <- !is.na(pos)
      if (any(take)) v[take] <- d[[cc]][pos[take]]
      svals[[cc]][[k]] <- v
    }
  }
  ord <- order(rep(seq_len(nu), times = ng), rep(seq_len(ng), each = nu))
  label <- themes$themes$label[match(t, themes$themes$theme_id)]
  out <- list(
    theme_id = rep(t, nu * ng),
    label = rep(label, nu * ng),
    gene = rep(genes, times = ng)[ord],
    group = rep(seq_len(ng), each = nu)[ord],
    evidence_members = unlist(em, use.names = FALSE)[ord],
    n_evidence_members = unlist(nem, use.names = FALSE)[ord],
    group_status = rep(status, each = nu)[ord]
  )
  for (cc in scols) out[[cc]] <- do.call(c, svals[[cc]])[ord]
  out$stat_status <- as.character(unlist(sstat, use.names = FALSE))[ord]
  out
}
