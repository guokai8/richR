# richSensitivity() (SPEC-RICHR-005 REQ-010 .. REQ-014): an ORA grid over
# user-named backgrounds x maxGSSize, one independent enrich() run per cell,
# states through bindRich(), a descriptive stability summary and its design.
# Fixture F-GRID (inline, seed 1) is the generator of ledger L-M5.

fgrid <- function() {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(1)
  U <- sprintf("G%04d", 1:2000)
  an <- data.frame(
    GeneID = c(unlist(lapply(1:6, function(i) sample(U, 600))),
               unlist(lapply(1:20, function(i) sample(U, 50)))),
    Term = c(rep(sprintf("B%d", 1:6), each = 600), rep(sprintf("S%02d", 1:20), each = 50)))
  q <- unique(c(an$GeneID[an$Term == "S01"][1:12], sample(U, 100)))
  list(U = U, an = an, q = q, bk = list(annotation = NULL, measured = U[1:1800]))
}

grid_cells <- function(bk, bounds = c(500, 1000)) {
  cells <- expand.grid(bound = bounds, background = names(bk), stringsAsFactors = FALSE)
  cells$cell <- paste0(cells$background, "/", format(cells$bound, scientific = FALSE, trim = TRUE))
  cells
}

# The state bindRich() assigns, computed here from one direct run's records.
direct_state <- function(r, id) {
  if (id %in% r@result$Annot) return("reported")
  if (id %in% testedTerms(r)$Annot) return("tested_not_reported")
  if (enrichmentInfo(r)$testing$zero_overlap_count > 0) "untested_or_zero_overlap" else "not_tested"
}

test_that("AC-010: every cell is an independent enrich() run on its own setting", {
  # Catches: a grid that filters one run instead of re-running (L-M5), cells
  # that do not vary only the stated axis, and per-cell states that disagree
  # with the engine's tested record.
  f <- fgrid()
  g <- richSensitivity(f$q, f$an, backgrounds = f$bk, maxGSSize = c(500, 1000))
  expect_s3_class(g, "richSensitivity")
  cells <- c("annotation/500", "annotation/1000", "measured/500", "measured/1000")
  expect_identical(g$cells$cell, cells)
  expect_identical(names(g$results), cells)
  expect_identical(g$cells$background_size, c(NA, NA, 1800L, 1800L))
  gc <- grid_cells(f$bk)
  for (i in seq_len(nrow(gc))) {
    d <- enrich(f$q, f$an, universe = f$bk[[gc$background[i]]], maxGSSize = gc$bound[i])
    r <- g$results[[gc$cell[i]]]
    expect_identical(r@result, d@result, info = gc$cell[i])
    expect_identical(testedTerms(r), testedTerms(d), info = gc$cell[i])
    expect_identical(enrichmentInfo(r)$testing, enrichmentInfo(d)$testing, info = gc$cell[i])
    t <- enrichmentInfo(d)$testing
    expect_identical(g$cells$effective_N[i], t$effective_N)
    expect_identical(g$cells$effective_n[i], t$effective_n)
    expect_identical(g$cells$tested_count[i], t$tested_count)
  }
  st <- g$states
  expect_identical(names(st)[1:4], c("group", "background", "maxGSSize", "Annot"))
  b500 <- st[st$group == "annotation/500" & grepl("^B", st$Annot), ]
  b1000 <- st[st$group == "annotation/1000" & grepl("^B", st$Annot), ]
  expect_identical(nrow(b500), 6L)
  expect_true(all(b500$state == "tested_not_reported"))
  expect_true(all(b1000$state %in% c("reported", "tested_not_reported")))
  t500 <- testedTerms(g$results[["annotation/500"]])
  t1000 <- testedTerms(g$results[["annotation/1000"]])
  expect_identical(t500["S01", "Padj"], t1000["S01", "Padj"])
})

test_that("AC-011: invalid calls error early; failed cells are kept and warnings recorded", {
  # Catches: a failed setting silently dropped, one bad level aborting the
  # grid, an all-failed grid returned as success, warnings swallowed or
  # repeated per cell, ambiguous cell ids, a measured-N request with no
  # measured set, a GSEA ranking run as if it were a gene list, and numeric
  # identifiers refused as if they were a ranking.
  f <- fgrid()
  rs <- function(...) richSensitivity(f$q, f$an, ...)
  expect_error(rs(backgrounds = unname(f$bk)), "`backgrounds`")
  expect_error(rs(backgrounds = list(a = NULL, a = f$U)), "`backgrounds`.*\"a\"")
  expect_error(rs(backgrounds = stats::setNames(list(NULL), " ")), "`backgrounds`")
  expect_error(rs(backgrounds = list(`a/b` = NULL)), "`backgrounds`.*a/b")
  expect_error(rs(backgrounds = list(a = 1:10)), "`backgrounds`.*\"a\"")
  for (bad in list(c(500, 500), 0, 500.5, numeric(0))) {
    expect_error(rs(backgrounds = f$bk, maxGSSize = bad), "`maxGSSize`")
  }
  expect_error(rs(backgrounds = f$bk, universe = f$U), "`universe`")
  expect_error(rs(backgrounds = f$bk, maxGSSize = c(500, 1000), maxGSSize = 10), "maxGSSize")
  expect_error(rs(backgrounds = f$bk, universe_policy = "measured"),
               "universe_policy.*annotation")
  expect_error(richSensitivity(c(A = 1, B = 2), f$an, backgrounds = f$bk), "`x`.*ranking")
  # Unnamed numbers are identifiers (e.g. integer Entrez IDs), not a ranking.
  ids_err <- tryCatch(richSensitivity(1:5, f$an, backgrounds = f$bk), error = conditionMessage)
  expect_match(ids_err, "`x`.*identifiers.*character")
  expect_no_match(ids_err, "ranking")

  bk3 <- c(f$bk, list(none = c("NOPE1", "NOPE2")))
  g <- rs(backgrounds = bk3)
  expect_identical(nrow(g$cells), 6L)
  none <- g$cells$background == "none"
  expect_identical(g$cells$execution[none], c("failed", "failed"))
  expect_match(g$cells$error[none], "universe has no annotated genes")
  expect_true(all(is.na(g$cells$effective_N[none])))
  expect_true(all(vapply(g$results[g$cells$cell[none]], is.null, logical(1))))
  expect_true(all(g$states$execution[g$states$background == "none"] == "failed"))
  expect_true(all(is.na(g$states$state[g$states$background == "none"])))
  g4 <- rs(backgrounds = f$bk)
  expect_identical(g$results[!none], g4$results)
  expect_identical(g$cells$execution[!none], g4$cells$execution)

  seen <- character()
  gd <- withCallingHandlers(
    richSensitivity(c(f$q, f$q[1]), f$an, backgrounds = f$bk),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_identical(sum(grepl("Duplicated genes", seen)), 1L)
  expect_true(all(vapply(gd$cells$warnings, function(w) any(grepl("Duplicated genes", w)),
                         logical(1))))

  expect_error(rs(backgrounds = list(none = c("NOPE1", "NOPE2"))),
               "none/500.*universe has no annotated genes.*none/1000")
})

test_that("AC-012: the summary counts every state from the cells' own records", {
  # Catches: a stability fraction over the wrong denominator, untested settings
  # counted as tested and absent, a possibly-tested zero-overlap term counted
  # as untested (or the reverse), failed cells counted as ran, identifiers no
  # cell tested listed with zero counts or zero-filled ranges, and score ranges
  # taken from anything but the cells' records.
  f <- fgrid()
  check_summary <- function(g, bk, ...) {
    gc <- grid_cells(bk)
    direct <- lapply(seq_len(nrow(gc)), function(i) tryCatch(
      enrich(f$q, f$an, universe = bk[[gc$background[i]]], maxGSSize = gc$bound[i], ...),
      error = function(e) NULL))
    ok <- !vapply(direct, is.null, logical(1))
    s <- g$summary
    expect_identical(names(s), c("Annot", "Term", "cells", "ran", "tested", "reported",
                                 "not_tested", "undetermined", "support", "padj_min",
                                 "padj_max", "fold_min", "fold_max", "rank_best",
                                 "rank_worst"))
    for (k in seq_len(nrow(s))) {
      id <- s$Annot[k]
      states <- vapply(direct[ok], direct_state, character(1), id = id)
      tested_cells <- direct[ok][states %in% c("reported", "tested_not_reported")]
      vals <- function(col) unlist(lapply(tested_cells, function(r) {
        t <- testedTerms(r)
        t[[col]][t$Annot == id]
      }))
      ranks <- unlist(lapply(tested_cells, function(r) {
        t <- testedTerms(r)
        rank(t$Pvalue, ties.method = "min")[t$Annot == id]
      }))
      rng <- function(v, fun) if (length(v)) fun(v) else NA_real_
      expect_identical(s$cells[k], nrow(gc), info = id)
      expect_identical(s$ran[k], sum(ok), info = id)
      expect_identical(s$tested[k], length(tested_cells), info = id)
      expect_identical(s$reported[k], sum(states == "reported"), info = id)
      expect_identical(s$not_tested[k], sum(states == "not_tested"), info = id)
      expect_identical(s$undetermined[k], sum(states == "untested_or_zero_overlap"), info = id)
      expect_identical(s$padj_min[k], rng(vals("Padj"), min), info = id)
      expect_identical(s$padj_max[k], rng(vals("Padj"), max), info = id)
      expect_identical(s$fold_min[k], rng(vals("FoldEnrichment"), min), info = id)
      expect_identical(s$fold_max[k], rng(vals("FoldEnrichment"), max), info = id)
      expect_identical(s$rank_best[k], as.integer(rng(ranks, min)), info = id)
      expect_identical(s$rank_worst[k], as.integer(rng(ranks, max)), info = id)
      want_support <- if (s$tested[k] == 0L) NA_real_ else s$reported[k] / s$tested[k]
      expect_identical(s$support[k], want_support, info = id)
    }
    s
  }

  bk3 <- c(f$bk, list(none = c("NOPE1", "NOPE2")))
  s <- check_summary(richSensitivity(f$q, f$an, backgrounds = bk3), bk3)
  expect_true(all(s$cells == 6L) && all(s$ran == 4L))
  s01 <- s[s$Annot == "S01", ]
  expect_identical(c(s01$tested, s01$reported, s01$support), c(4, 4, 1))
  b1 <- s[s$Annot == "B1", ]
  expect_identical(c(b1$tested, b1$not_tested), c(4L, 0L))
  expect_identical(b1$support, b1$reported / 4)
  zero <- s[s$tested > 0 & s$reported == 0, ]
  expect_gt(nrow(zero), 0L)
  expect_true(all(zero$support == 0))
  # An identifier no cell tests is absent from the summary, not listed with
  # zero counts or zero-filled ranges: with maxGSSize = 500 alone no cell tests
  # the 600-gene B sets, and the summary holds exactly the tested identifiers.
  g5 <- richSensitivity(f$q, f$an, backgrounds = f$bk, maxGSSize = 500)
  tested5 <- unique(unlist(lapply(g5$results, function(r) testedTerms(r)$Annot)))
  expect_true(any(grepl("^B", tested5)))
  expect_setequal(g5$summary$Annot, tested5)
  expect_true(all(g5$summary$tested >= 1L))
  expect_false(anyNA(g5$summary[c("padj_min", "padj_max", "fold_min", "fold_max",
                                  "rank_best", "rank_worst", "support")]))

  # Under upstream policy size bounds act after adjustment. Positive-overlap
  # large terms remain tested, with unchanged adjusted p-values across bounds.
  expect_true(all(vapply(g5$results, function(r)
    enrichmentInfo(r)$testing$zero_overlap_count == 0L, logical(1))))
})

test_that("AC-013: the design states what was varied and reproduces the grid", {
  # Catches: a sensitivity report that cannot state what was varied, held
  # fixed or counted, a grid that cannot be reproduced from its own record,
  # and a combined statistic slipped into the summary.
  f <- fgrid()
  g <- richSensitivity(f$q, f$an, backgrounds = f$bk, maxGSSize = c(500, 1000), pvalue = 0.01)
  d <- g$design
  expect_identical(d$backgrounds, f$bk)
  expect_identical(d$background_sizes, c(annotation = NA_integer_, measured = 1800L))
  expect_identical(d$maxGSSize, c(500, 1000))
  expect_identical(d$engine, "enrich")
  expect_identical(d$arguments, list(pvalue = 0.01))
  expect_identical(d$criterion$state, "reported")
  expect_identical(d$criterion[c("pvalue", "padj", "padj.method", "family_policy",
                                 "universe_policy")],
                   list(pvalue = 0.01, padj = NULL, padj.method = "BH",
                        family_policy = "positive_overlap_before_reporting_filters", universe_policy = "annotated"))
  expect_true(all(c("ran", "tested", "support") %in% names(d$definitions)))
  expect_match(d$interpretation, "descriptive")
  expect_match(d$interpretation, "not a p-value")
  expect_identical(d$richR_version, as.character(utils::packageVersion("richR")))
  expect_false(any(grepl("pvalue|p_value|combined|mean", names(g$summary), ignore.case = TRUE)))

  again <- do.call(richSensitivity, c(list(f$q, f$an, d$backgrounds, d$maxGSSize), d$arguments))
  expect_identical(again$cells, g$cells)
  expect_identical(again$states, g$states)
  expect_identical(again$summary, g$summary)
})

test_that("a measured-background grid over Annot preserves explicitly requested settings", {
  f <- fgrid()
  an <- methods::new("Annot", species="x", anntype="custom", keytype="SYMBOL", annot=f$an)
  bk <- list(measured=f$U)
  ga <- richSensitivity(f$q, an, backgrounds=bk, universe_policy="measured")
  gd <- richSensitivity(f$q, f$an, backgrounds=bk, universe_policy="measured")
  expect_identical(ga$design$arguments$universe_policy,"measured")
  expect_identical(ga$design$criterion$universe_policy,"measured")
  expect_identical(ga$cells$tested_count,gd$cells$tested_count)
  for(cell in ga$cells$cell) expect_identical(ga$results[[cell]]@tested,gd$results[[cell]]@tested)
  expect_identical(ga$summary,gd$summary)
})
