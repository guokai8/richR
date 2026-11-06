# Analysis audit (SPEC-RICHR-005 REQ-006 .. REQ-009): the engines keep the
# mapping record of their input and the per-identifier input counts, and
# mappingAudit() reads one decision per identifier from the recorded metadata.
# Fixture F-AUD is inline; F-ORA-C, F-GSEA and M4's legacy fixtures are files.

faud_annot <- function() {
  data.frame(GeneID = paste0("G", 1:20), Term = rep(paste0("T", 1:4), each = 5))
}
faud_qa <- c("G1", "G2", "G3", "G1", "G6", NA, "", "X9", "G19")
faud_ua <- c(paste0("G", 1:18), "U1", "U1", NA)
faud_map <- function() {
  data.frame(from = c("g1", "g2", "g3a", "g3a", "x9", "g6", "g6b", "g19"),
             to = c("G1", "G2", "G3", "G4", "X9", "G6", "G6", "G19"))
}
faud_mapped_query <- c("g1", "g2", "g3a", "x9", "zz", "g6", "g6b", "g19")
faud_enrich <- function(q, universe = faud_ua, annot = faud_annot(), ...) {
  suppressWarnings(enrich(q, annot, pvalue = 1, padj = NULL, minSize = 1, minGSSize = 1,
                          universe = universe, ...))
}
info_input <- function(r) enrichmentInfo(r)$input

# readRich() refuses any attribute other than names (M4, L-R1).
no_class_or_row_names <- function(x) {
  if (any(c("class", "row.names") %in% names(attributes(x)))) return(FALSE)
  if (is.list(x)) return(all(vapply(x, no_class_or_row_names, logical(1))))
  TRUE
}

identity_mapped <- function(x) {
  ids <- if (is.numeric(x)) names(x) else x
  mapGenes(x, data.frame(from = ids, to = ids), "from", "to")
}

gsea_fixture_ranking <- function() {
  rk <- read.csv(test_path("fixtures", "gsea_ranked_a.csv"), stringsAsFactors = FALSE)
  stats::setNames(rk$score, rk$GeneID)
}
gsea_fixture_annot <- function() {
  read.csv(test_path("fixtures", "gsea_annot.csv"), colClasses = "character")
}
run_gsea <- function(s) {
  suppressWarnings(richGSEA(s, gsea_fixture_annot(), minSize = 3, maxSize = 20,
                            seed = 1, nproc = 1))
}

test_that("AC-005: engines keep a matching mapping record and the input occurrences", {
  # Catches: provenance of an upstream conversion lost at an engine (including
  # richKEGG), a stale record attached to a modified input, the record changing
  # statistics, duplicated identifiers unrecoverable (L-M3), and a record lost,
  # retyped or stored in a shape readRich() refuses (L-R1).
  v <- mapGenes(faud_mapped_query, faud_map(), "from", "to")
  rec <- attr(v, "richR_mapping")
  r <- faud_enrich(v)
  expect_identical(info_input(r)$mapping, rec)
  expect_true(no_class_or_row_names(info_input(r)$mapping))

  annot_obj <- new("Annot", species = "x", anntype = "custom", keytype = "SYMBOL",
                   annot = faud_annot())
  expect_identical(info_input(faud_enrich(v, annot = annot_obj))$mapping, rec)

  oq <- read.csv(test_path("fixtures", "ora_query.csv"), colClasses = "character")$GeneID
  ka <- read.csv(test_path("fixtures", "kegg_annot.csv"), colClasses = "character")
  mq <- identity_mapped(oq)
  rk <- suppressMessages(richKEGG(mq, ka, pvalue = 1, padj = NULL, minSize = 1,
                                  minGSSize = 1))
  expect_identical(info_input(rk)$mapping, attr(mq, "richR_mapping"))

  plain <- faud_enrich(as.vector(v))
  expect_identical(r@result, plain@result)
  expect_null(info_input(plain)$mapping)
  expect_no_warning(u <- enrich(unique(v), faud_annot(), pvalue = 1, padj = NULL,
                                minSize = 1, minGSSize = 1, universe = faud_ua))
  expect_null(info_input(u)$mapping)
  seen <- character()
  rq <- withCallingHandlers(
    enrich(faud_qa, faud_annot(), pvalue = 1, padj = NULL, minSize = 1, minGSSize = 1,
           universe = faud_ua),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("mapping", seen)))
  expect_null(info_input(rq)$mapping)
  expect_identical(info_input(rq)$ids, c("G1", "G2", "G3", "G6", "X9", "G19"))
  expect_identical(info_input(rq)$occurrences, c(2L, 1L, 1L, 1L, 1L, 1L))

  v2 <- v
  v2[1] <- "G5"
  expect_warning(r2 <- enrich(v2, faud_annot(), pvalue = 1, padj = NULL, minSize = 1,
                              minGSSize = 1, universe = faud_ua),
                 "mapping record does not match the submitted input")
  expect_identical(info_input(r2)$mapping,
                   list(available = FALSE,
                        reason = "mapping record does not match the submitted input"))

  p <- tempfile("m5-audit-")
  writeRich(r, p)
  y <- readRich(p)
  expect_identical(enrichmentInfo(y), enrichmentInfo(r))
  expect_identical(enrichmentInfo(y)$input$mapping$decisions$candidates[[5]], character(0))

  skip_if_not_installed("fgsea")
  s <- gsea_fixture_ranking()
  ms <- identity_mapped(s)
  g <- run_gsea(ms)
  expect_identical(info_input(g)$mapping, attr(ms, "richR_mapping"))
  expect_true(no_class_or_row_names(info_input(g)$mapping))
  pg <- tempfile("m5-audit-gsea-")
  writeRich(g, pg)
  expect_identical(enrichmentInfo(readRich(pg)), enrichmentInfo(g))
})

test_that("AC-006: the plain ORA audit gives every identifier one reasoned decision", {
  # Catches: an excluded identifier with no visible reason, an annotated gene
  # outside the universe reported as unannotated (or the reverse), duplicated
  # or invalid entries missing, counts that do not reconcile with n and N.
  r <- faud_enrich(faud_qa)
  a <- mappingAudit(r)
  expect_identical(names(a), c("role", "input_id", "mapping", "id", "occurrences",
                               "decision", "reason"))
  q <- a[a$role == "query", ]
  b <- a[a$role == "background", ]
  ids <- q[!is.na(q$id), ]
  expect_identical(ids$id, c("G1", "G2", "G3", "G6", "X9", "G19"))
  expect_identical(ids$decision, c(rep("effective", 4), "not_annotated", "outside_background"))
  expect_identical(ids$occurrences, c(2L, 1L, 1L, 1L, 1L, 1L))
  expect_identical(ids$input_id, ids$id)
  expect_true(all(is.na(q$mapping)))
  expect_match(ids$reason[ids$id == "X9"], "outside the supplied universe")
  inv <- q[is.na(q$id), ]
  expect_identical(nrow(inv), 2L)
  expect_true(all(is.na(inv$input_id)))
  expect_identical(inv$decision, c("invalid", "invalid"))
  expect_identical(inv$occurrences, c(1L, 1L))
  expect_identical(nrow(b), 19L)
  expect_setequal(b$id[b$decision == "effective"], paste0("G", 1:18))
  expect_identical(b$id[b$decision == "not_annotated"], "U1")

  info <- enrichmentInfo(r)
  expect_identical(length(unique(q$id[q$decision == "effective"])), info$testing$effective_n)
  expect_identical(info$testing$effective_n, 4L)
  expect_identical(sum(b$decision == "effective"), info$testing$effective_N)
  expect_identical(info$testing$effective_N, 18L)
  expect_identical(sum(q$occurrences, na.rm = TRUE), info$input$original_count)
  expect_identical(info$input$original_count, 9L)

  n <- mappingAudit(faud_enrich(faud_qa, universe = NULL))
  nq <- n[n$role == "query" & !is.na(n$id), ]
  expect_identical(nq$decision[nq$id == "G19"], "effective")
  expect_identical(nq$decision[nq$id == "X9"], "not_annotated")
  expect_false(grepl("universe", nq$reason[nq$id == "X9"]))
  nb <- n[n$role == "background", ]
  expect_identical(nrow(nb), 20L)
  expect_true(all(nb$decision == "effective"))
})

test_that("AC-007: the joined ORA audit connects every original id to the engine decision", {
  # Catches: identifiers lost at conversion invisible in the result's audit, a
  # merged entry double-counted against n or the input count, and an audit
  # that cannot connect a user's original identifier to the engine decision.
  v <- mapGenes(faud_mapped_query, faud_map(), "from", "to")
  r <- faud_enrich(v)
  a <- mappingAudit(r)
  q <- a[a$role == "query", ]
  expect_identical(q$input_id, faud_mapped_query)
  expect_identical(q$mapping, c("mapped", "mapped", "ambiguous", "mapped", "unmapped",
                                "mapped", "merged", "mapped"))
  expect_identical(q$id, c("G1", "G2", NA, "X9", NA, "G6", "G6", "G19"))
  expect_identical(q$decision, c("effective", "effective", "not_submitted", "not_annotated",
                                 "not_submitted", "effective", "effective",
                                 "outside_background"))
  expect_identical(q$occurrences, c(1L, 1L, NA, 1L, NA, 1L, NA, 1L))
  info <- enrichmentInfo(r)
  expect_identical(sum(q$occurrences, na.rm = TRUE), info$input$original_count)
  expect_identical(info$input$original_count, 5L)
  expect_identical(length(unique(q$id[q$decision == "effective"])), info$testing$effective_n)
  expect_identical(info$testing$effective_n, 3L)
  expect_identical(sum(a$role == "background" & a$decision == "effective"),
                   info$testing$effective_N)
  expect_identical(mappingAudit(v), list2DF(attr(v, "richR_mapping")$decisions))
})

test_that("AC-008: the GSEA audit lists every ranked gene, its score and ties", {
  # Catches: unannotated ranked genes reported as excluded (they stay in the
  # ranking), tied genes left uninspectable, scores misattributed to genes, and
  # conversion losses missing from a GSEA audit.
  skip_if_not_installed("fgsea")
  s <- gsea_fixture_ranking()
  s[c("G10", "G11", "G12")] <- s[["G10"]]
  s <- c(s, Z1 = 0.0123, Z2 = -0.0456)
  g <- run_gsea(s)
  a <- mappingAudit(g)
  info <- enrichmentInfo(g)
  expect_identical(names(a), c("role", "input_id", "mapping", "id", "occurrences",
                               "decision", "reason", "score", "tied"))
  expect_true(all(a$role == "ranked"))
  expect_identical(a$id, info$input$ids)
  expect_identical(a$occurrences, rep(1L, length(s)))
  annotated <- unique(gsea_fixture_annot()$GeneID)
  expect_identical(a$decision, ifelse(a$id %in% annotated, "ranked_annotated",
                                      "ranked_not_annotated"))
  expect_identical(a$decision[a$id %in% c("Z1", "Z2")], rep("ranked_not_annotated", 2))
  expect_identical(a$score, unname(s[a$id]))
  tie <- duplicated(s) | duplicated(s, fromLast = TRUE)
  expect_identical(a$tied, unname(tie[match(a$id, names(s))]))
  expect_identical(sum(a$tied), 3L)
  expect_identical(length(unique(a$id[a$decision == "ranked_annotated"])),
                   info$input$mapped_count)

  base <- gsea_fixture_ranking()
  low <- stats::setNames(base, tolower(names(base)))
  low <- c(low, extra1 = base[["G03"]], extra2 = base[["G07"]] + 0.5)
  tab <- data.frame(from = c(tolower(names(base)), "extra1", "extra2"),
                    to = c(names(base), "G03", "G07"))
  ms <- mapGenes(low, tab, "from", "to")
  gm <- run_gsea(ms)
  j <- mappingAudit(gm)
  dec <- attr(ms, "richR_mapping")$decisions
  expect_identical(j$input_id, names(low))
  expect_identical(j$mapping, dec$decision)
  merged <- which(j$mapping == "merged")
  expect_identical(j$input_id[merged], "extra1")
  expect_identical(j$id[merged], "G03")
  expect_identical(j$decision[merged], j$decision[j$input_id == "g03"])
  expect_identical(j$occurrences[merged], NA_integer_)
  conflict <- which(j$mapping == "conflict")
  expect_identical(j$input_id[conflict], c("g07", "extra2"))
  expect_identical(j$decision[conflict], rep("not_submitted", 2))
  expect_true(all(is.na(j$id[conflict])))
  expect_identical(j$score[conflict], unname(low[c("g07", "extra2")]))
  expect_true(all(is.na(j$tied[conflict])))
  expect_identical(j$occurrences[j$mapping == "mapped"], rep(1L, sum(j$mapping == "mapped")))
  expect_identical(length(unique(j$id[j$decision == "ranked_annotated"])),
                   enrichmentInfo(gm)$input$mapped_count)
})

test_that("AC-009: legacy metadata gives a plain audit; missing metadata gives none", {
  # Catches: an audit reconstructed from @gene or @result (fabricated
  # provenance), an error on pre-0.2.0 objects or on a malformed mapping
  # record, and a missing occurrence count filled with 1.
  lr <- readRDS(test_path("fixtures", "legacy_richResult_fca5276.rds"))
  mr <- enrichmentInfo(lr)
  expect_no_warning(a <- mappingAudit(lr))
  q <- a[a$role == "query" & !is.na(a$id), ]
  expect_identical(q$id, mr$input$ids)
  expect_identical(q$id[q$decision == "effective"], mr$input$effective_ids)
  expect_true(all(is.na(q$occurrences)))
  expect_true(all(is.na(a$mapping)))

  lg <- readRDS(test_path("fixtures", "legacy_GSEAResult_fca5276.rds"))
  mg <- enrichmentInfo(lg)
  expect_no_warning(ag <- mappingAudit(lg))
  expect_identical(ag$id, mg$input$ids)
  expect_identical(ag$id[ag$decision == "ranked_annotated"], mg$input$mapped_ids)
  expect_true(all(is.na(ag$occurrences)))
  expect_true(all(is.na(ag$mapping)))

  r <- faud_enrich(faud_qa)
  coerced <- methods::as(r@result, "richResult")
  stripped <- r
  attr(stripped, "richR_metadata") <- NULL
  for (obj in list(coerced, stripped, c("G1", "G2"))) {
    expect_warning(out <- mappingAudit(obj), "mapping audit unavailable: .+")
    expect_null(out)
  }
  # A vector whose richR_mapping attribute is not a mapGenes() record gets the
  # same unavailable result, not an internal error.
  for (rec in list(list(foo = 1), "x", list(decisions = list(decision = 1)))) {
    expect_warning(out <- mappingAudit(structure("a", richR_mapping = rec)),
                   "mapping audit unavailable: the mapping record is malformed")
    expect_null(out)
  }
})
