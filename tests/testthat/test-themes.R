# richThemes(): theme grouping of reported terms (SPEC-RICHR-007 M6.2a).
# Fixtures and their role checks live in helper-themes.R.

# ---- AC-001: grouping set, kappa similarity, average linkage ----------------

test_that("every reported F-TH term lands in exactly one kappa theme (catches terms dropped or themes filtered as in richCluster(), and a similarity other than kappa on evidence genes)", {
  .th_fth_roles()
  x <- .th_fth()
  th <- richThemes(list(g1 = x$g1, g2 = x$g2), cutoff = 0.5)

  expect_identical(.th_partition(th), c("ABC", "DE", "F", "G"))
  m <- unique(th$members[, c("theme_id", "Annot")])
  expect_identical(sort(m$Annot, method = "radix"), LETTERS[1:7])
  expect_identical(th$settings$background_size, 27L)

  s <- .th_fth_sets(); U <- unique(unlist(s))
  sim <- th$similarity
  expect_identical(nrow(sim), 4L)
  pair <- paste(sim$identifier_a, sim$identifier_b)
  expect_setequal(pair, c("A B", "A C", "B C", "D E"))
  own <- mapply(function(a, b) .th_kappa(s[[a]], s[[b]], U), sim$identifier_a, sim$identifier_b)
  expect_lt(max(abs(sim$kappa - own)), 1e-12)
  expect_identical(sim$n_shared[match(c("A B", "A C", "B C", "D E"), pair)], c(7L, 7L, 6L, 7L))

  th9 <- richThemes(list(g1 = x$g1, g2 = x$g2), cutoff = 0.9)
  expect_identical(.th_partition(th9), LETTERS[1:7])
  expect_identical(nrow(th9$themes), 7L)
})

test_that("F-LINK is cut with average linkage, not single or complete (catches a different linkage changing which terms share a theme)", {
  .th_flink_roles()
  th <- richThemes(.th_flink(), cutoff = 0.5)
  expect_identical(.th_partition(th), c("PQR", "T", "Z"))
})

test_that("the cut is inclusive: terms whose kappa equals cutoff share a theme (catches richCluster()'s strict kappa > cutoff, which splits them)", {
  u <- function(i) sprintf("u%02d", i)
  # P-Q kappa is exactly 0.5 on the background u01..u04; R shares nothing
  expect_identical(.th_kappa(u(1), u(1:2), u(1:4)), 0.5)
  pqr <- .th_res(c("P", "Q", "R"), list(u(1), u(1:2), u(3:4)), c(0.01, 0.02, 0.03),
                 c("p", "q", "r"))
  th <- richThemes(pqr, cutoff = 0.5)
  expect_identical(.th_partition(th), c("PQ", "R"))
  expect_identical(th$similarity$kappa, 0.5)

  # K and M have identical similarity genes (kappa 1): they join at cutoff = 1
  x <- .th_fth2()
  expect_identical(.th_partition(richThemes(list(q1 = x$q1, q2 = x$q2), cutoff = 1)),
                   c("KM", "L", "N"))
})

test_that("evidence genes are split on the result's own sep (catches a non-comma sep read as one gene per term, which changes every theme)", {
  G <- sprintf("G%02d", 1:80)
  q <- G[c(1:6, 11:14, 17)]
  comma <- .th_enrich(q)
  slash <- .th_enrich(q, sep = "/")
  expect_identical(slash@sep, "/")
  expect_identical(slash@result$GeneID[slash@result$Annot == "K"], "G01/G02/G03/G04/G05/G06")
  a <- richThemes(list(q1 = comma))
  s <- richThemes(list(q1 = slash))
  expect_identical(s$genes$gene[s$genes$Annot == "K"], sprintf("G%02d", 1:6))
  expect_identical(s$themes, a$themes)
  expect_identical(s$genes, a$genes)
  expect_identical(s$similarity, a$similarity)
})

test_that("terms tested but reported nowhere are not grouped (catches tested-only terms shown as if they were findings)", {
  x <- .th_fth2()
  .th_fth2_roles(x$q1, x$q2)
  th <- richThemes(list(q1 = x$q1, q2 = x$q2))
  expect_false("W" %in% th$members$Annot)
  expect_identical(.th_partition(th), c("KM", "L", "N"))
  expect_identical(th$settings$background_size, 16L)
})

# ---- AC-002: invariance to group and row order ------------------------------

test_that("themes, ids, representatives and labels do not change with group or row order (catches a label or theme id taken from the input order, which breaks saved links to a theme)", {
  x <- .th_fth()
  rev_rows <- function(o) { o@result <- o@result[rev(seq_len(nrow(o@result))), ]; o }
  a <- richThemes(list(g1 = x$g1, g2 = x$g2))
  b <- richThemes(list(g2 = x$g2, g1 = x$g1))
  c <- richThemes(list(g1 = rev_rows(x$g1), g2 = rev_rows(x$g2)))
  # the D/E label really is at stake: g2 holds another Term for D
  expect_identical(x$g2@result$Term[x$g2@result$Annot == "D"], "death receptor (g2 label)")

  expect_identical(a$themes, b$themes)
  expect_identical(a$themes, c$themes)
  expect_identical(a$similarity, b$similarity)
  expect_identical(a$similarity, c$similarity)

  norm_m <- function(m) {
    m$group <- as.character(m$group)
    m <- m[order(m$theme_id, m$Annot, m$group, method = "radix"), ]
    rownames(m) <- NULL
    m
  }
  norm_g <- function(g) {
    g$group <- as.character(g$group)
    g <- g[order(g$theme_id, g$Annot, g$group, g$gene, method = "radix"), ]
    rownames(g) <- NULL
    g
  }
  expect_identical(norm_m(a$members), norm_m(b$members))
  expect_identical(norm_m(a$members), norm_m(c$members))
  expect_identical(norm_g(a$genes), norm_g(b$genes))
  expect_identical(norm_g(a$genes), norm_g(c$genes))
})

# ---- AC-003: the themes table ------------------------------------------------

test_that("the representative is the most central member, not the smallest p-value, and no theme carries a p-value (catches a representative's p-value presented as the theme's)", {
  x <- .th_fth()
  th <- richThemes(list(g1 = x$g1, g2 = x$g2))
  tt <- th$themes
  expect_identical(names(tt), c("theme_id", "label", "representative", "representative_term",
                                "n_terms", "n_groups_reported", "n_similarity_genes",
                                "mean_within_kappa"))
  expect_identical(tt$theme_id, sprintf("theme_%03d", 1:4))
  expect_identical(tt$representative, c("A", "D", "F", "G"))
  expect_identical(tt$n_terms, c(3L, 2L, 1L, 1L))

  # A is the representative although it has the largest Padj of A, B, C
  padj <- stats::setNames(x$g1@result$Padj, x$g1@result$Annot)
  expect_true(padj[["A"]] > padj[["B"]] && padj[["A"]] > padj[["C"]])
  s <- .th_fth_sets(); U <- unique(unlist(s))
  own <- c(.th_kappa(s$A, s$B, U), .th_kappa(s$A, s$C, U), .th_kappa(s$B, s$C, U))
  expect_lt(abs(tt$mean_within_kappa[1] - mean(own)), 1e-12)
  expect_identical(tt$mean_within_kappa[3:4], c(NA_real_, NA_real_))

  expect_identical(tt$label[1], "alpha")
  # D's Term in the radix-first reporting group g1, although g2 holds another
  expect_identical(tt$label[2], "death receptor")
  expect_identical(tt$representative_term[4], NA_character_)
  expect_identical(tt$label[4], "G")
  expect_identical(tt$n_groups_reported, c(1L, 2L, 1L, 1L))
  expect_identical(tt$n_similarity_genes, c(10L, 9L, 4L, 4L))

  pcols <- c("Pvalue", "Padj", "pval", "padj", "NES", "EnrichmentScore")
  for (el in setdiff(names(th), c("members", "settings"))) {
    expect_length(intersect(names(th[[el]]), pcols), 0L)
  }
  expect_length(intersect(names(th$settings), pcols), 0L)
  pv <- c(x$g1@result$Pvalue, x$g1@result$Padj, x$g2@result$Pvalue, x$g2@result$Padj)
  nums <- unlist(tt[vapply(tt, is.numeric, logical(1))])
  expect_false(any(nums[!is.na(nums)] %in% pv))

  th9 <- richThemes(list(g1 = x$g1, g2 = x$g2), cutoff = 0.9)
  expect_identical(utils::tail(th9$themes$theme_id, 1), "theme_007")
})

test_that("themes are numbered by size before representative (catches theme ids assigned in representative order, which moves every theme_id link when a larger theme has a later representative)", {
  y <- function(i) sprintf("y%02d", i)
  a <- sprintf("a%02d", 1:3)
  # Y and Z join (kappa 0.58); the singleton A sorts before Y
  expect_gte(.th_kappa(y(1:6), y(c(1:5, 7)), c(a, y(1:7))), 0.5)
  ayz <- .th_res(c("A", "Y", "Z"), list(a, y(1:6), y(c(1:5, 7))), c(0.01, 0.02, 0.03),
                 c("a term", "y term", "z term"))
  tt <- richThemes(ayz)$themes
  expect_identical(tt$theme_id, c("theme_001", "theme_002"))
  expect_identical(tt$representative, c("Y", "A"))
  expect_identical(tt$n_terms, c(2L, 1L))
})

test_that("a blank Term is not a label: the next reporting group's Term, else the identifier, is used (catches a whitespace or empty label shown for a theme)", {
  y <- function(i) sprintf("y%02d", i)
  g1 <- .th_res(c("A", "B"), list(y(1:4), y(11:14)), c(0.01, 0.02), c("  ", ""))
  g2 <- .th_res("A", list(y(1:4)), 0.03, "alpha")
  tt <- richThemes(list(g1 = g1, g2 = g2))$themes
  expect_identical(tt$representative, c("A", "B"))
  # A: g1 is radix-first but its Term is blank, so g2's Term is used
  expect_identical(tt$representative_term, c("alpha", NA))
  expect_identical(tt$label, c("alpha", "B"))
})

# ---- AC-004: reversibility to bindRich() rows and the groups' records -------

test_that("members are the bindRich() rows and genes come from each group's own record (catches a theme that cannot be traced to its rows, and genes of a reporting group shown for a tested-not-reported member)", {
  x <- .th_fth2()
  .th_fth2_roles(x$q1, x$q2)
  b <- bindRich(list(q1 = x$q1, q2 = x$q2))
  th <- richThemes(list(q1 = x$q1, q2 = x$q2))

  expect_identical(names(th$members), c("theme_id", "is_representative", names(b)))
  mm <- th$members[, names(b)]
  bb <- b[b$Annot %in% c("K", "M", "L", "N"), ]
  mm <- mm[order(mm$Annot, mm$group, method = "radix"), ]
  bb <- bb[order(bb$Annot, bb$group, method = "radix"), ]
  rownames(mm) <- NULL; rownames(bb) <- NULL
  expect_identical(mm, bb)

  g <- th$genes
  expect_identical(names(g), c("theme_id", "group", "Annot", "gene", "evidence_type", "state"))
  k1 <- g[g$Annot == "K" & g$group == "q1", ]
  k2 <- g[g$Annot == "K" & g$group == "q2", ]
  expect_identical(k1$gene, sprintf("G%02d", 1:6))
  expect_true(all(k1$state == "reported"))
  expect_identical(k2$gene, c("G07", "G08"))
  expect_true(all(k2$state == "tested_not_reported"))
  expect_identical(nrow(g[g$Annot %in% c("M", "L") & g$group == "q2", ]), 0L)
  expect_identical(nrow(g[g$Annot == "N" & g$group == "q1", ]), 0L)
  expect_false("W" %in% g$Annot)
  expect_true(all(g$evidence_type == "ora_overlap"))

  km <- th$similarity
  expect_identical(nrow(km), 1L)
  expect_identical(c(km$identifier_a, km$identifier_b), c("K", "M"))
  expect_identical(km$kappa, 1)
  expect_identical(km$n_shared, 6L)

  expect_setequal(unique(th$members$Annot[th$members$is_representative]), c("K", "L", "N"))
  expect_false(any(th$members$is_representative[th$members$Annot == "M"]))

  created <- list(themes = names(th$themes), genes = names(th$genes),
                  similarity = names(th$similarity), ranking = names(th$ranking),
                  members = c("theme_id", "is_representative"))
  for (el in c("themes", "members", "genes", "similarity", "ranking")) {
    tab <- th[[el]]
    cls <- vapply(tab, richR:::.rich_io_column_class, character(1))
    expect_false(anyNA(cls), info = paste(el, paste(names(cls)[is.na(cls)], collapse = ",")))
    own <- vapply(tab[created[[el]]], function(v) class(v)[1], character(1))
    expect_true(all(own %in% c("character", "integer", "numeric", "logical", "factor")),
                info = el)
  }
})

test_that("GSEA leading-edge genes are read the same from both storage forms and the ranked input is kept (catches leading edges lost or mis-split, genes invented for an empty edge, and a dropped or re-sorted ranking)", {
  skip_if_not_installed("fgsea")
  t1 <- .th_fgsea(TRUE); t0 <- .th_fgsea(FALSE)
  for (o in c(t1, t0)) expect_true(.th_t5_empty(o))
  expect_type(t0$a@result$leadingEdge, "list")
  expect_type(t1$a@result$leadingEdge, "character")

  th1 <- richThemes(list(a = t1$a, b = t1$b))
  th0 <- richThemes(list(a = t0$a, b = t0$b))
  expect_gt(nrow(th1$genes), 0L)
  expect_identical(th1$genes, th0$genes)
  expect_true(all(th1$genes$evidence_type == "gsea_leading_edge"))
  expect_false("T5" %in% th1$genes$pathway)
  expect_identical(th1$themes$n_terms[th1$themes$representative == "T5"], 1L)
  expect_identical(sum(th1$members$pathway == "T5" & th1$members$theme_id ==
                         th1$themes$theme_id[th1$themes$representative == "T5"]), 2L)

  for (gname in c("a", "b")) {
    r <- th1$ranking[th1$ranking$group == gname, ]
    expect_identical(stats::setNames(r$rank_score, r$gene), t1[[gname]]@input)
  }
  expect_identical(names(th1$ranking), c("group", "gene", "rank_score"))
  expect_type(th0$members$leadingEdge, "list")
})

test_that("the ranking keeps the stored order of an input that is not sorted (catches a ranking re-sorted in either direction)", {
  skip_if_not_installed("fgsea")
  t1 <- .th_fgsea(TRUE)
  a <- t1$a
  n <- length(a@input)
  a@input <- a@input[c(2L, 1L, seq.int(3L, n))]
  expect_false(identical(a@input, sort(a@input)))
  expect_false(identical(a@input, sort(a@input, decreasing = TRUE)))
  rk <- richThemes(list(a = a, b = t1$b))$ranking
  r <- rk[rk$group == "a", ]
  expect_identical(stats::setNames(r$rank_score, r$gene), a@input)
})

# ---- AC-006: settings and a failed group ------------------------------------

test_that("settings record the method, threshold and each group's own source run, and a failed group is kept as failed (catches unreadable provenance, one group's record stored under another, and a NULL group dropped)", {
  x <- .th_fth2()
  .th_fth2_roles(x$q1, x$q2)
  th <- richThemes(list(q1 = x$q1, q2 = x$q2, f = NULL))
  th0 <- richThemes(list(q1 = x$q1, q2 = x$q2))
  s <- th$settings
  expect_identical(s$method, "kappa")
  expect_identical(s$distance, "1 - kappa")
  expect_identical(s$linkage, "average")
  expect_identical(s$cutoff, 0.5)
  expect_identical(s$groups, c("q1", "q2", "f"))
  expect_identical(s$background_size, 16L)
  expect_identical(s$sources$q1, enrichmentInfo(x$q1))
  expect_identical(s$sources$q2, enrichmentInfo(x$q2))
  expect_false(identical(s$sources$q1, s$sources$q2))
  expect_identical(s$sources$f, "failed")
  expect_identical(names(s$sources), c("q1", "q2", "f"))
  expect_identical(s$richR_version, as.character(utils::packageVersion("richR")))
  for (f in c("grouped_on", "similarity_genes", "background", "representative", "label")) {
    expect_true(is.character(s[[f]]) && length(s[[f]]) == 1L && nzchar(s[[f]]), info = f)
  }

  fr <- th$members[th$members$group == "f", ]
  expect_gt(nrow(fr), 0L)
  expect_true(all(fr$execution == "failed"))
  expect_true(all(is.na(fr$state)))
  expect_identical(levels(th$members$group), c("q1", "q2", "f"))
  expect_identical(th$themes, th0$themes)
})

# ---- AC-009: edge inputs ----------------------------------------------------

test_that("a run reporting nothing gives empty typed tables, not an error (catches an empty analysis turned into an error or into untyped columns SPEC-RICHR-008 cannot store)", {
  z <- .th_fzero()
  bz <- bindRich(list(z = z))
  expect_identical(nrow(z@result), 0L)
  expect_identical(bz$state, "tested_not_reported")

  types <- list(
    themes = c(theme_id = "character", label = "character", representative = "character",
               representative_term = "character", n_terms = "integer",
               n_groups_reported = "integer", n_similarity_genes = "integer",
               mean_within_kappa = "numeric"),
    genes = c(theme_id = "character", group = "factor", Annot = "character",
              gene = "character", evidence_type = "character", state = "character"),
    similarity = c(theme_id = "character", identifier_a = "character",
                   identifier_b = "character", kappa = "numeric", n_shared = "integer"),
    ranking = c(group = "factor", gene = "character", rank_score = "numeric"))
  for (input in list(list(z = z, groups = "Group 1", obj = z),
                     list(z = z, groups = c("z", "f"), obj = list(z = z, f = NULL)))) {
    th <- richThemes(input$obj)
    expect_s3_class(th, "richThemes")
    expect_identical(names(th), c("themes", "members", "genes", "similarity", "ranking", "settings"))
    for (el in names(types)) {
      expect_identical(nrow(th[[el]]), 0L)
      expect_identical(vapply(th[[el]], function(v) class(v)[1], ""), types[[el]], info = el)
    }
    expect_identical(nrow(th$members), 0L)
    expect_identical(names(th$members),
                     c("theme_id", "is_representative",
                       names(bindRich(if (is.list(input$obj)) input$obj else list(input$obj)))))
    expect_identical(class(th$members$theme_id), "character")
    expect_identical(class(th$members$is_representative), "logical")
    expect_identical(th$settings$groups, input$groups)
    expect_identical(th$settings$background_size, 0L)
    expect_identical(levels(th$genes$group), input$groups)
  }
  expect_identical(richThemes(.th_fth()$g1)$settings$groups, "Group 1")
})

test_that("empty evidence makes a singleton without genes and identical evidence keeps the documented kappa 0 (catches genes fabricated for an empty list and the degenerate case silently repaired)", {
  x <- .th_fth()
  s <- .th_fth_sets()
  g1e <- .th_res(c("A", "B", "C", "D", "E", "F", "H"),
                 c(s[c("A", "B", "C", "D", "E", "F")], list(character(0))),
                 c(0.04, 0.001, 0.01, 0.02, 0.03, 0.005, 0.01),
                 c("alpha", "beta", "gamma", "death receptor", "epsilon", "phi", "eta"))
  expect_identical(g1e@result$GeneID[7], "")
  th <- richThemes(g1e)
  expect_identical(th$themes$n_terms[th$themes$representative == "H"], 1L)
  expect_false("H" %in% th$genes$Annot)
  expect_identical(.th_partition(th), c("ABC", "DE", "F", "H"))

  dg <- .th_res(c("X", "Y"), list(c("h01", "h02", "h03"), c("h01", "h02", "h03")),
                c(0.01, 0.02), c("x", "y"))
  td <- richThemes(dg)
  expect_identical(.th_partition(td), c("X", "Y"))
  expect_identical(td$settings$background_size, 3L)
})

test_that("a legacy group without a tested table is grouped from its result and its absent terms stay unavailable (catches invented states or genes for a legacy group)", {
  x <- .th_fth2()
  leg <- .th_fleg()
  th <- richThemes(list(q1 = x$q1, leg = leg))
  lm <- th$members[th$members$group == "leg", ]
  st <- stats::setNames(lm$state, lm$Annot)
  expect_identical(unname(st[c("K", "P", "M", "L")]),
                   c("reported", "reported", "unavailable", "unavailable"))
  expect_true(all(c("K", "P") %in% th$members$Annot))
  lg <- th$genes[th$genes$group == "leg", ]
  expect_identical(paste(lg$Annot, lg$gene), c("K G01", "K G02", "K G03", "P G50", "P G51"))
  expect_true(all(lg$state == "reported"))
})

test_that("mixed ORA/GSEA lists and mixed leadingEdge forms stop with bindRich()'s error (catches ORA and GSEA terms mixed in one theme set)", {
  skip_if_not_installed("fgsea")
  q1 <- .th_fth2()$q1
  t1 <- .th_fgsea(TRUE); t0 <- .th_fgsea(FALSE)
  for (bad in list(list(q1, t1$a), list(a = t0$a, b = t1$b))) {
    msg <- tryCatch(bindRich(bad), error = conditionMessage)
    expect_true(is.character(msg) && nzchar(msg))
    expect_error(richThemes(bad), msg, fixed = TRUE)
  }
})

test_that("an invalid cutoff is refused and cutoff = 1 is accepted (catches a threshold outside (0, 1] silently accepted)", {
  g1 <- .th_fth()$g1
  for (cu in list(0, 1.5, NA, Inf, c(0.3, 0.4), "0.5", NA_real_, -0.2)) {
    expect_error(richThemes(g1, cutoff = cu), "cutoff", fixed = TRUE)
  }
  th <- richThemes(g1, cutoff = 1)
  expect_identical(th$settings$cutoff, 1)
  expect_identical(nrow(th$themes), 6L)
})
