# bindRich() (SPEC-RICHR-004 REQ-008 .. REQ-010): one long table with an
# explicit execution state and inferential state per term per group. Fixtures
# F-ORA and F-GSEA are built in helper-tested.R with the public API.

# Named state vector of one group, keyed by identifier.
bind_states <- function(b, group, id = "Annot") {
  rows <- b[b$group == group, , drop = FALSE]
  stats::setNames(rows$state, rows[[id]])
}

# Every row takes Term and its statistics only from its own group's record:
# @result for "reported", the tested table for "tested_not_reported", and NA
# (never 0) in every other state and on every row of a failed group.
check_bound_values <- function(b, objs, id = "Annot") {
  value_cols <- setdiff(names(b), c("group", id, "state", "execution"))
  for (g in names(objs)) {
    rows <- b[b$group == g, , drop = FALSE]
    x <- objs[[g]]
    for (i in seq_len(nrow(rows))) {
      state <- rows$state[[i]]
      src <- if (identical(state, "reported")) {
        as.data.frame(x@result)
      } else if (identical(state, "tested_not_reported")) {
        testedTerms(x)
      }
      if (!is.null(src) && identical(id, "pathway") && !"Term" %in% names(src)) {
        src$Term <- as.character(src$pathway)
      }
      for (col in value_cols) {
        got <- rows[[col]][[i]]
        if (is.null(src) || !col %in% names(src)) {
          expect_true(all(is.na(got)), info = paste(g, rows[[id]][[i]], col))
        } else {
          want <- src[[col]][[match(rows[[id]][[i]], src[[id]])]]
          expect_identical(got, want, info = paste(g, rows[[id]][[i]], col))
        }
      }
    }
  }
}

# ---- States (AC-007, AC-005 bindRich leg) ------------------------------------

test_that("bindRich() states per group what was reported, tested, untested or unknowable", {
  # Catches missing or untested values written as 0, execution and inferential
  # states mixed in one column, a failed group dropped, an empty run confused
  # with a failed one, "not tested" asserted where the record cannot know it,
  # and rows fabricated for identifiers no group recorded (T4).
  full <- .fora_reference()
  .fora_roles(testedTerms(full))
  # Explicit historical fixture: this engine retains its positive-overlap family.
  all <- full
  m <- enrichmentInfo(all)
  m$testing$zero_overlap_count <- 1L
  m$testing$tested_count <- as.integer(nrow(all@tested) + 1L)
  m$testing$family_policy <- "explicit_historical_fixture"
  attr(all, "richR_metadata") <- m
  expect_gt(enrichmentInfo(all)$testing$zero_overlap_count, 0L)
  legacy <- as(full@result, "richResult")
  other <- suppressMessages(enrich(c("G18", "G19"), .fora_annot(), pvalue = 1,
                                   padj = NULL, minSize = 2, minGSSize = 3))
  expect_identical(other@result$Annot, "T3")
  expect_identical(testedTerms(other)$Annot, "T3")
  none <- .fora_reference(pvalue = 1e-12)
  expect_identical(nrow(none@result), 0L)
  nometa <- full
  attr(nometa, "richR_metadata") <- NULL

  objs <- list(full = full, all = all, legacy = legacy, other = other,
               none = none, nometa = nometa, failed = NULL)
  b <- bindRich(objs)

  expect_identical(class(b), "data.frame")
  expect_identical(nrow(b), 7L * 5L)
  expect_setequal(unique(b$Annot), c("T1", "T2", "T3", "T4", "T5"))
  expect_identical(levels(b$group), names(objs))
  execution <- vapply(names(objs),
                      function(g) unique(b$execution[b$group == g]), "")
  expect_identical(execution, c(full = "success", all = "success",
                                legacy = "success", other = "success",
                                none = "empty", nometa = "success",
                                failed = "failed"))

  ids <- c("T1", "T2", "T3", "T5")
  tnr <- "tested_not_reported"
  expected <- list(
    full = c("reported", tnr, "not_tested", tnr),
    all = c("reported", tnr, "untested_or_zero_overlap", tnr),
    legacy = c("reported", "unavailable", "unavailable", "unavailable"),
    other = c("not_tested", "not_tested", "reported", "not_tested"),
    none = c(tnr, tnr, "not_tested", tnr),
    nometa = c("reported", tnr, "unavailable", tnr),
    failed = rep(NA_character_, 4)
  )
  for (g in names(expected)) {
    expect_identical(unname(bind_states(b, g)[ids]), expected[[g]], info = g)
  }
  check_bound_values(b, objs[names(objs) != "failed"])
  failed <- b[b$group == "failed",
              setdiff(names(b), c("group", "Annot", "execution"))]
  expect_true(all(is.na(failed)))
})

test_that("bindRich() has a fixed column contract across ORA entry points", {
  # Catches columns that differ between implementations, a statistic filled
  # with 0 for a group that never had it, and the tested table's reported
  # flag leaking into the bound table.
  e <- .fora_reference()
  k <- suppressMessages(richKEGG(
    .char_read("ora_query.csv")$GeneID,
    .char_read("kegg_annot.csv", character = TRUE),
    pvalue = 0.05, padj = NULL, minSize = 2, minGSSize = 3
  ))
  l <- as(e@result, "richResult")
  b <- bindRich(list(e = e, k = k, l = l))
  expect_identical(names(b), c(
    "group", "Annot", "Term", "state", "execution", "Annotated",
    "Significant", "RichFactor", "FoldEnrichment", "zscore", "Pvalue", "Padj",
    "GeneID", "ko", "Level3", "Level2", "Level1"
  ))
  el <- b$group %in% c("e", "l")
  for (col in c("ko", "Level3", "Level2", "Level1")) {
    expect_identical(b[[col]][el], rep(NA_character_, sum(el)), info = col)
  }
  expect_identical(b$RichFactor[b$group == "l"],
                   rep(NA_real_, sum(b$group == "l")))
  check_bound_values(b, list(e = e, k = k, l = l))
})

test_that("bindRich() binds derived views by their current @result", {
  # Catches a filtered view whose dropped terms are reported as untested, a
  # grouped @result bound differently from its plain rows, and a view without
  # its label column bound as if complete.
  r <- .fora_reference()
  f <- bindRich(list(f = dplyr::filter(r, Annot != "T1")))
  expect_identical(bind_states(f, "f")[["T1"]], "tested_not_reported")
  expect_identical(f$execution[f$Annot == "T1"], "empty")

  grouped <- r
  grouped@result <- dplyr::group_by(r@result, Term)
  expect_identical(bindRich(list(g = grouped)),
                   bindRich(list(g = r)))
  expect_error(bindRich(list(s = dplyr::select(r, -Term))), "'s'.*Term")
})

# ---- GSEA and labels (AC-008) -------------------------------------------------

test_that("bindRich() binds GSEA results by pathway and never zero-fills untested pathways", {
  # Catches a GSEA path that binds by label, GSEA "not tested" collapsed into
  # zero, and a legacy GSEA record reported as complete.
  skip_if_not_installed("fgsea")
  a <- .fgsea_run(minSize = 3, maxSize = 20)
  bb <- .fgsea_run("gsea_ranked_b.csv", minSize = 6, maxSize = 20)
  expect_true("T4" %in% testedTerms(a)$pathway)
  expect_false("T4" %in% testedTerms(bb)$pathway)
  legacy <- a
  attr(legacy, "tested") <- NULL
  attr(legacy, "richR_metadata") <- NULL

  objs <- list(a = a, b = bb, legacy = legacy)
  out <- bindRich(objs)
  expect_identical(names(out)[1:5],
                   c("group", "pathway", "Term", "state", "execution"))
  t4b <- out[out$group == "b" & out$pathway == "T4", ]
  expect_identical(t4b$state, "not_tested")
  expect_true(all(is.na(t4b[, c("NES", "pval", "padj", "Term")])))

  leg <- bind_states(out, "legacy", id = "pathway")
  absent <- setdiff(names(leg), legacy@result$pathway)
  expect_gt(length(absent), 0L)
  expect_true(all(leg[absent] == "unavailable"))
  expect_true(all(leg[legacy@result$pathway] == "reported"))
  check_bound_values(out, objs, id = "pathway")
})

test_that("bindRich() keeps identifiers whose label is blank or missing", {
  # Catches rows dropped or keyed by label when a Term is "" or NA.
  r <- .fora_reference()
  x <- as(data.frame(Annot = c("X1", "X2"), Term = c("", NA),
                     Annotated = c(4L, 4L), Significant = c(2L, 2L),
                     Pvalue = c(0.01, 0.02), Padj = c(0.02, 0.04),
                     GeneID = c("G01,G02", "G03,G04")), "richResult")
  out <- bindRich(list(r = r, x = x))
  xr <- out[out$group == "x" & out$Annot %in% c("X1", "X2"), ]
  expect_identical(xr$Annot, c("X1", "X2"))
  expect_identical(xr$Term, c("", NA))
  expect_identical(xr$state, c("reported", "reported"))
  expect_identical(bind_states(out, "r")[c("X1", "X2")],
                   c(X1 = "not_tested", X2 = "not_tested"))
})

test_that("bindRich() fills list and factor columns from each group's own record", {
  # Catches leading edges (a list column, table = FALSE) or factor values
  # shifted between rows or groups, filled for untested rows, or factor levels
  # of the later groups dropped.
  skip_if_not_installed("fgsea")
  a <- .fgsea_run(minSize = 3, maxSize = 20, table = FALSE)
  bb <- .fgsea_run("gsea_ranked_b.csv", minSize = 6, maxSize = 20, table = FALSE)
  expect_type(a@result$leadingEdge, "list")
  objs <- list(a = a, b = bb)
  out <- bindRich(objs)
  expect_type(out$leadingEdge, "list")
  expect_identical(out$leadingEdge[out$group == "b" & out$pathway == "T4"],
                   list(NA_character_))
  check_bound_values(out, objs, id = "pathway")

  r1 <- .fora_reference()
  r2 <- .fora_reference(pvalue = 1)
  r1@result$kind <- factor("x", levels = c("x", "y"))
  r1@tested$kind <- factor(rep("y", nrow(r1@tested)), levels = c("x", "y"))
  r2@result$kind <- factor(rep("z", nrow(r2@result)), levels = "z")
  out <- bindRich(list(a = r1, b = r2))
  expect_s3_class(out$kind, "factor")
  expect_identical(levels(out$kind), c("x", "y", "z"))
  st <- paste(out$group, out$state)
  want <- c("a reported" = "x", "a tested_not_reported" = "y",
            "b reported" = "z")[st]
  expect_identical(as.character(out$kind), unname(want))
})

# ---- Input validation (AC-009) ------------------------------------------------

test_that("bindRich() preserves legacy default group labels and rejects ambiguous names", {
  # Catches group labels or name errors that disagree with compareResult()
  # for the same list.
  r1 <- .fora_reference()
  expect_identical(levels(bindRich(list(r1, r1))$group),
                   c("Group 1", "Group 2"))
  for (bad in list(list(a = r1, a = r1), list(a = r1, " " = r1))) {
    expect_error(bindRich(bad), "duplicate|blank")
  }
})

test_that("bindRich() rejects inputs it cannot bind without guessing", {
  # Catches a data.frame bound as if it carried a record, silent coercion of
  # conflicting columns, a view without its identifier bound as if complete,
  # and an all-failed batch returned as an empty success.
  r1 <- .fora_reference()
  expect_error(bindRich(list(a = r1, b = r1@result)), "'b'")
  expect_error(bindRich(list(a = r1, b = dplyr::mutate(r1, GeneID = 1))),
               "GeneID.*'a'.*'b'")
  expect_error(bindRich(list(a = dplyr::rename(r1, ID = Annot))),
               "'a'.*Annot")
  expect_error(bindRich(list(a = NULL, b = NULL)), "NULL")
  expect_error(bindRich(list()), "non-empty list")

  # Integer and double are one numeric type: widened, not rejected.
  w <- bindRich(list(a = r1,
                     b = dplyr::mutate(r1, Annotated = as.double(Annotated))))
  expect_type(w$Annotated, "double")
  expect_identical(w$Annotated[w$group == "a" & w$Annot == "T1"],
                   as.double(r1@result$Annotated))
})

test_that("bindRich() rejects a list mixing ORA and GSEA results", {
  # Catches ORA and GSEA scores silently placed in one table (proposal 9.4).
  skip_if_not_installed("fgsea")
  r1 <- .fora_reference()
  g1 <- .fgsea_run(minSize = 3, maxSize = 20)
  expect_error(bindRich(list(a = r1, b = g1)), "ORA.*'a'.*GSEA.*'b'")
})
