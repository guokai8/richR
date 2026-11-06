# richReport() study report (SPEC-RICHR-008 AC-001 ... AC-008). Fixtures and
# the base-R HTML readers are in helper-report.R; each test asserts its
# fixture roles first and names the defect it catches.

test_that("report tables, saved results and checksums reproduce what the report shows (AC-001: lost levels, list columns, a failed group shifting indices)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  .rp_fq_roles(q$q1, q$q2)
  p <- file.path(base, "p")
  m <- .rp_main(p, q = q)
  q1 <- q$q1; q2 <- q$q2; x <- list(q1 = q1, q2 = q2, f = NULL)

  ref <- file.path(base, "ref")
  writeRich(q1, ref)
  result_files <- sort(list.files(ref), method = "radix")
  tables <- c("comparison", "audit_group_1", "audit_group_2", "themes", "theme_members",
              "theme_genes", "theme_similarity", "theme_ranking", "theme_gene_view",
              "theme_stats_group_1", "sensitivity_cells", "sensitivity_summary",
              "sensitivity_states")
  expected <- sort(c("report.html", "manifest.json", "rerun.R",
                     file.path("results/group_1", result_files),
                     file.path("results/group_2", result_files),
                     paste0("tables/", tables, ".tsv")), method = "radix")
  expect_identical(.rp_files(p), expected)

  expect_identical(readRich(file.path(p, "results/group_1")), q1)
  expect_identical(readRich(file.path(p, "results/group_2")), q2)

  b <- bindRich(x)
  expect_identical(levels(b$group), c("q1", "q2", "f"))
  expect_identical(.rp_decode(p, "comparison"), b)
  expect_identical(.rp_decode(p, "audit_group_1"), suppressWarnings(mappingAudit(q1)))
  expect_identical(.rp_decode(p, "audit_group_2"), suppressWarnings(mappingAudit(q2)))
  expect_identical(.rp_decode(p, "themes")$label[1], "Edited theme")
  expect_identical(.rp_decode(p, "themes"), m$th$themes)
  expect_identical(.rp_decode(p, "theme_members"), m$th$members)
  expect_identical(.rp_decode(p, "theme_genes"), m$th$genes)
  expect_identical(.rp_decode(p, "theme_similarity"), m$th$similarity)
  expect_identical(.rp_decode(p, "theme_ranking"), m$th$ranking)
  expect_identical(.rp_decode(p, "theme_gene_view"), themeGenes(m$th, stats = .rp_st()))
  expect_identical(.rp_decode(p, "theme_stats_group_1"), .rp_st()$q1)
  expect_true(is.list(m$ss$cells$warnings))
  expect_identical(.rp_decode(p, "sensitivity_cells"), m$ss$cells)
  expect_identical(.rp_decode(p, "sensitivity_summary"), m$ss$summary)
  expect_identical(.rp_decode(p, "sensitivity_states"), m$ss$states)

  mf <- .rp_manifest(p)
  expect_identical(mf$format, "richR-report")
  expect_identical(vapply(mf$groups, `[[`, "", "name"), c("q1", "q2", "f"))
  expect_identical(lapply(mf$groups, `[[`, "folder"),
                   list("results/group_1", "results/group_2", NULL))
  listed <- names(mf$checksums)
  expect_setequal(listed, setdiff(expected, c("report.html", "manifest.json")))
  for (f in listed) {
    expect_identical(digest::digest(file = file.path(p, f), algo = "sha256"),
                     mf$checksums[[f]], label = f)
  }
})

test_that("GSEA report tables keep the leadingEdge list column and every audit (AC-001, AC-003: GSEA evidence lost in the tables)", {
  skip_if_not_installed("fgsea")
  base <- .rp_tmp()
  gs <- .th_fgsea(FALSE)
  expect_true(is.list(gs$a@result$leadingEdge))
  pg <- file.path(base, "pg")
  richReport(gs, pg, themes = richThemes(gs))
  cmp <- .rp_decode(pg, "comparison")
  expect_identical(cmp, bindRich(gs))
  expect_true(is.list(cmp$leadingEdge))
  expect_identical(.rp_decode(pg, "audit_group_1"), suppressWarnings(mappingAudit(gs$a)))
  expect_identical(.rp_decode(pg, "audit_group_2"), suppressWarnings(mappingAudit(gs$b)))
  expect_identical(readRich(file.path(pg, "results/group_2")), gs$b)
})

test_that("report statements lead back to their term, run and annotation (AC-002: unanchored rows, numbers in absent states, a theme shown as a test, a nested record dumped or blank)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  .rp_fq_roles(q$q1, q$q2)
  q1 <- q$q1; q2 <- q$q2
  leg <- .rp_fleg(q1)
  bu <- bindRich(list(q2 = q2, leg = leg))
  expect_identical(bu$state[bu$Annot == "N" & bu$group == "leg"], "unavailable")
  expect_identical(bu$state[bu$Annot == "M" & bu$group == "q2"], "not_tested")

  p <- file.path(base, "p")
  .rp_main(p, q = q)
  h <- .rp_html(p)
  ids <- .rp_ids(h)
  expect_true(all(.rp_anchor(1, q1@result$Annot) %in% ids))
  expect_true(all(.rp_anchor(2, q2@result$Annot) %in% ids))
  local <- sub("^#", "", grep("^#", .rp_hrefs(h), value = TRUE))
  expect_gt(length(local), 0)
  expect_true(all(local %in% ids))

  # Comparison cells against an expectation built from bindRich().
  expect_cells <- function(html, x, groups) {
    b <- bindRich(x)
    rows <- .rp_table_rows(.rp_section(html, "comparison"))
    expect_identical(rows[[1]], c("identifier", groups))
    body <- rows[-1]
    rep_ids <- sort(unique(b$Annot[b$state %in% "reported"]), method = "radix")
    expect_identical(vapply(body, `[`, "", 1), rep_ids)
    for (r in body) {
      for (k in seq_along(groups)) {
        s <- b[b$Annot == r[1] & as.character(b$group) == groups[k], ]
        want <- if (is.null(x[[groups[k]]])) "failed"
                else if (s$state %in% c("reported", "tested_not_reported"))
                  sprintf("%s (%s)", s$state, sprintf("%.4g", s$Padj))
                else s$state
        expect_identical(r[k + 1], want, label = paste(r[1], groups[k]))
      }
    }
    body
  }
  body <- expect_cells(h, list(q1 = q1, q2 = q2, f = NULL), c("q1", "q2", "f"))
  cell <- function(body, id) body[[which(vapply(body, `[`, "", 1) == id)]][-1]
  expect_identical(cell(body, "K"), c("reported (0.0001042)", "tested_not_reported (0.1967)", "failed"))
  expect_identical(cell(body, "N"), c("not_tested", "reported (6.823e-06)", "failed"))
  expect_true(.rp_has(h, "4 identifiers reported in at least one group, of 5 in tables/comparison.tsv."))
  for (s in c("A term significant in one group and not in another is not evidence of a difference between the groups.",
              "Theme grouping summarizes the listed terms; it is not a statistical test, and no theme-level p-value exists.",
              "Evidence genes are ORA overlap genes or GSEA leading-edge genes; they are not shown to be causal drivers.",
              "Related or overlapping gene sets are not independent replications.")) {
    expect_true(.rp_has(h, s), label = s)
  }

  # Edited label shown, members linked to their q1 rows.
  h3 <- .rp_match_all(h, "<h3 id=\"theme_001\">.*?</h3>")
  expect_identical(trimws(.rp_text(h3)), "Edited theme (theme_001)")
  t1 <- .rp_section(h, "theme_001", until = "<h3")
  members <- unique(.rp_th(q1, q2)$members$Annot[.rp_th(q1, q2)$members$theme_id == "theme_001"])
  expect_setequal(sub("^#", "", .rp_hrefs(t1)), .rp_anchor(1, members))
  # A missing fold change (q2 has no statistics) is shown as NA, never 0.
  gv <- .rp_table_rows(t1)
  hit <- Filter(function(r) identical(r[match(c("gene", "group"), gv[[1]])], c("G01", "q2")),
                gv[-1])
  expect_length(hit, 1L)
  expect_identical(hit[[1]][match(c("log2FC", "stat_status"), gv[[1]])],
                   c("NA", "no_statistics_for_group"))

  # Provenance and the record summary rules.
  for (o in list(q1, q2)) {
    expect_true(grepl(enrichmentInfo(o)$provenance$annotation$fingerprint$membership, h,
                      fixed = TRUE))
  }
  g1 <- .rp_section(h, "group-1")
  expect_identical(.rp_record_value(g1, "testing$effective_n"), "11")
  expect_identical(.rp_record_value(g1, "background$effective_ids"), "character vector of 76 values")
  expect_identical(.rp_record_value(g1, "input$unmapped_ids"), "character vector of 0 values")
  expect_identical(.rp_record_value(g1, "options$padj"), "NULL")
  # Every reported row with every column, in result order.
  rt <- .rp_table_rows(g1, "Reported terms: ")
  expect_identical(rt[[1]], names(q1@result))
  expect_identical(vapply(rt[-1], `[`, "", match("Annot", rt[[1]])), q1@result$Annot)

  pl <- file.path(base, "pl")
  richReport(list(q1 = q1, leg = leg), pl)
  expect_true(.rp_has(.rp_section(.rp_html(pl), "group-2"),
                      "Analysis record unavailable: metadata not recorded."))

  pu <- file.path(base, "pu")
  richReport(list(q2 = q2, leg = leg), pu)
  bu_rows <- expect_cells(.rp_html(pu), list(q2 = q2, leg = leg), c("q2", "leg"))
  expect_identical(cell(bu_rows, "N"), c("reported (6.823e-06)", "unavailable"))
  expect_identical(cell(bu_rows, "M"), c("not_tested", "reported (5.555e-05)"))
})

test_that("the mapping audit is shown as recorded, never fabricated or recomputed (AC-003: mapping evidence dropped, an audit invented for a legacy object, a mismatched record shown as complete)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  qm <- .rp_fqm()
  am <- suppressWarnings(mappingAudit(qm))
  expect_identical(nrow(am), 89L)
  expect_identical(sum(am$role == "query" & am$decision == "not_submitted"), 1L)

  p <- file.path(base, "p")
  richReport(list(qm = qm, q2 = q$q2), p)
  expect_identical(.rp_decode(p, "audit_group_1"), am)
  h <- .rp_html(p)
  g1 <- .rp_section(h, "group-1")
  expect_true(.rp_has(g1, "Mapping audit: 89 rows."))
  rows <- .rp_table_rows(g1, "Mapping audit: 89 rows.")[-1]
  got <- lapply(rows, function(r) r[1:3])
  tab <- table(am$role, am$decision)
  nz <- which(tab > 0, arr.ind = TRUE)
  want <- lapply(seq_len(nrow(nz)), function(i) {
    c(rownames(tab)[nz[i, 1]], colnames(tab)[nz[i, 2]], as.character(tab[nz[i, 1], nz[i, 2]]))
  })
  want <- want[order(vapply(want, `[`, "", 1), vapply(want, `[`, "", 2), method = "radix")]
  expect_identical(got, want)
  expect_identical(got, list(c("background", "effective", "76"), c("query", "effective", "12"),
                             c("query", "not_submitted", "1")))
  expect_true("tables/audit_group_1.tsv" %in% .rp_hrefs(g1))
  expect_true(.rp_has(g1, "Shown as recorded; not recomputed by richReport()."))
  expect_identical(.rp_record_value(g1, "input$mapping$decisions$candidates"), "list of 13 values")

  pl <- file.path(base, "pl")
  richReport(list(q1 = q$q1, leg = .rp_fleg(q$q1)), pl)
  expect_false(file.exists(file.path(pl, "tables/audit_group_2.tsv")))
  expect_true(.rp_has(.rp_section(.rp_html(pl), "group-2"),
                      "Mapping audit unavailable: metadata not recorded."))

  qx <- .rp_fqx()
  expect_warning(ax <- mappingAudit(qx), "mapping record is unavailable")
  expect_identical(nrow(ax), 87L)
  px <- file.path(base, "px")
  seen <- character(0)
  withCallingHandlers(richReport(list(qx = qx), px),
                      warning = function(w) seen <<- c(seen, conditionMessage(w)))
  expect_false(any(startsWith(seen, "mappingAudit()")))
  expect_identical(.rp_decode(px, "audit_group_1"), ax)
  hx <- .rp_html(px)
  expect_true(.rp_has(hx, "Mapping audit: 87 rows."))
  expect_true(.rp_has(hx, paste("Mapping audit note: mappingAudit(): the mapping record is unavailable",
                                "(mapping record does not match the submitted input); showing the",
                                "analysis input only.")))
})

test_that("GSEA audits are written per group (AC-003: GSEA mapping evidence dropped)", {
  skip_if_not_installed("fgsea")
  base <- .rp_tmp()
  gs <- .th_fgsea(FALSE)
  pg <- file.path(base, "pg")
  richReport(gs, pg)
  expect_identical(.rp_decode(pg, "audit_group_1"), suppressWarnings(mappingAudit(gs$a)))
  expect_true(.rp_has(.rp_html(pg), sprintf("Mapping audit: %d rows.",
                                            nrow(suppressWarnings(mappingAudit(gs$a))))))
})

test_that("figures draw every reported row up to figure.top and say what they draw (AC-004: re-thresholding at plot defaults, a misstated caption, a zero-row group crashing the report)", {
  skip_if_not(capabilities("png"), "the png device is unavailable")
  base <- .rp_tmp()
  fr <- .rp_ffig()
  expect_identical(nrow(fr@result), 12L)
  expect_identical(sum(fr@result$Pvalue > 0.05), 4L)
  expect_identical(fr@result$Annot, sprintf("S%02d", 1:12))

  p1 <- file.path(base, "p1"); richReport(fr, p1, figure.top = 30)
  p2 <- file.path(base, "p2"); richReport(fr, p2, figure.top = 5)
  expect_identical(.rp_captions(.rp_html(p1)), "12 of 12 reported terms; stored rows with valid drawing values; first 30 in result order; 0 omitted for missing, blank or non-finite drawing values; 0 truncated by top")
  expect_identical(.rp_captions(.rp_html(p2)), "5 of 12 reported terms; stored rows with valid drawing values; first 5 in result order; 0 omitted for missing, blank or non-finite drawing values; 7 truncated by top")
  f1 <- .rp_manifest(p1)$figures_drawn[[1]]
  f2 <- .rp_manifest(p2)$figures_drawn[[1]]
  expect_identical(c(f1$k, f1$n, f2$k, f2$n), c(12L, 12L, 5L, 12L))
  expect_identical(f2$rule, "stored rows with valid drawing values; first 5 in result order")
  expect_identical(f2[["function"]], "richDot")

  q <- .rp_fq()
  z <- .rp_fzero()
  expect_identical(nrow(z@result), 0L)
  pz <- file.path(base, "pz"); richReport(list(z = z, q1 = q$q1), pz)
  hz <- .rp_html(pz)
  expect_identical(length(.rp_match_all(hz, "data:image/png;base64,")), 1L)
  expect_identical(.rp_captions(hz)[1], "No figure: no reported term has a p-value.")

  p3 <- file.path(base, "p3"); richReport(fr, p3, figures = FALSE)
  h3 <- .rp_html(p3)
  expect_false(grepl("<img", h3, fixed = TRUE))
  expect_identical(.rp_captions(h3), "No figure: figures = FALSE.")

  pn <- file.path(base, "pn"); richReport(list(q1 = q$q1, f = NULL), pn)
  hn <- .rp_html(pn)
  expect_identical(.rp_captions(hn)[2], "No figure: the group failed.")
  expect_true(.rp_has(.rp_section(hn, "group-2"), "Execution: failed."))
})

test_that("GSEA figures follow richNES()'s order rule and caption it (AC-004: a GSEA caption that misstates the drawn rows)", {
  skip_if_not_installed("fgsea")
  skip_if_not(capabilities("png"), "the png device is unavailable")
  base <- .rp_tmp()
  gs <- .th_fgsea(FALSE)
  expect_identical(c(nrow(gs$a@result), nrow(gs$b@result)), c(5L, 5L))
  pg <- file.path(base, "pg")
  richReport(gs, pg, figure.top = 3)
  expect_identical(.rp_captions(.rp_html(pg)),
                   rep("3 of 5 reported terms; stored rows with finite NES; first 3 by NES, descending; 0 omitted for missing, blank or non-finite drawing values; 2 truncated by top", 2))
  expect_identical(.rp_manifest(pg)$figures_drawn[[2]][["function"]], "richNES")
})

test_that("sensitivity evidence is stored and shown as supplied (AC-005: stability evidence altered or dropped, NA shown as 0, universe vectors leaked, a non-plain design written)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  x <- list(q1 = q$q1, q2 = q$q2)
  ss <- .rp_fsens()
  expect_identical(nrow(ss$cells), 4L)
  expect_identical(ss$cells$background_size[ss$cells$cell == "all/10"], NA_integer_)

  p <- file.path(base, "p")
  richReport(x, p, sensitivity = ss)
  expect_identical(.rp_decode(p, "sensitivity_cells"), ss$cells)
  expect_identical(.rp_decode(p, "sensitivity_summary"), ss$summary)
  expect_identical(.rp_decode(p, "sensitivity_states"), ss$states)
  node <- .rp_manifest(p)$sensitivity_design
  d <- jsonlite::unserializeJSON(as.character(jsonlite::toJSON(
    node, auto_unbox = TRUE, null = "null", digits = I(17))))
  expect_identical(d$background_names, c("all", "part"))
  expect_identical(unname(d$background_sizes), c(NA, 60L))
  expect_identical(d$maxGSSize, c(10, 50))
  expect_identical(d$interpretation, ss$design$interpretation)
  mtext <- paste(readLines(file.path(p, "manifest.json"), warn = FALSE), collapse = "\n")
  expect_false(grepl("G(0[1-9]|[1-5][0-9]|60)\\b", mtext, perl = TRUE))

  h <- .rp_section(.rp_html(p), "sensitivity")
  expect_true(.rp_has(h, ss$design$interpretation))
  expect_true(.rp_has(h, "Cells: 4 rows."))
  expect_true(.rp_has(h, "Summary: 4 rows."))
  expect_true(.rp_has(h, "Shown as supplied; not recomputed by richReport()."))
  rows <- .rp_table_rows(h, "Cells: 4 rows.")
  col <- match("background_size", rows[[1]])
  row <- Filter(function(r) identical(r[1], "all/10"), rows[-1])[[1]]
  expect_identical(row[col], "NA")

  p2 <- file.path(base, "p2")
  richReport(x, p2)
  expect_false(any(grepl("sensitivity", .rp_files(p2))))
  expect_false(grepl("id=\"sensitivity\"", .rp_html(p2), fixed = TRUE))

  p3 <- file.path(base, "p3")
  expect_error(richReport(x, p3, sensitivity = ss$summary), "sensitivity")
  ss2 <- ss
  ss2$design$arguments$fun <- identity
  p4 <- file.path(base, "p4")
  expect_error(richReport(x, p4, sensitivity = ss2), "arguments")
  expect_false(file.exists(p3))
  expect_false(file.exists(p4))
})

test_that("report.html is inert and every hostile text is escaped (AC-006: script injection through titles, group names, terms, theme labels or genes; remote resources; broken links)", {
  base <- .rp_tmp()
  old <- Sys.getenv(c("http_proxy", "https_proxy"), unset = NA)
  on.exit({
    for (v in names(old)) if (is.na(old[[v]])) Sys.unsetenv(v) else do.call(Sys.setenv, as.list(old[v]))
  }, add = TRUE)
  Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9")

  H <- c("<script>alert(1)</script>", "a&b \"c\" 'd'", "<img src=x onerror=1>")
  q <- .rp_fq()
  qh <- q$q1
  qh@result$Term[1:3] <- H
  expect_identical(qh@result$Term[1:3], H)
  gene <- "<i>G01</i>"
  qk <- .rp_enrich(replace(.rp_q1_query(), 1, gene), annot = .rp_a80(c(G01 = gene)))
  expect_true(grepl(gene, qk@result$GeneID[qk@result$Annot == "K"], fixed = TRUE))
  thk <- richThemes(list(k = qk))
  thk$themes$label[1] <- "<script>x</script>"
  title <- "<script>t</script> & co"

  x <- stats::setNames(list(qh, q$q2), c("<b>g1</b>", "q2"))
  on_img <- file.path(base, "on"); off_img <- file.path(base, "off"); pk <- file.path(base, "pk")
  richReport(x, on_img, title = title, figures = TRUE)
  richReport(x, off_img, title = title, figures = FALSE)
  richReport(list(k = qk), pk, themes = thk)

  allowed <- c("html", "head", "meta", "title", "style", "body", "h1", "h2", "h3", "p",
               "table", "tr", "th", "td", "figure", "img", "figcaption", "a", "code")
  for (p in c(on_img, off_img, pk)) {
    h <- .rp_html(p)
    tags <- .rp_match_all(h, "<[A-Za-z][^>]*>")
    expect_gt(length(tags), 0)
    expect_true(all(grepl("^<[a-z0-9]+( [a-z-]+=\"[^\"<>]*\")*>$", tags)), label = p)
    expect_true(all(sub("^<([A-Za-z][A-Za-z0-9]*).*$", "\\1", tags) %in% allowed), label = p)
    expect_identical(length(.rp_match_all(h, "<img")), length(.rp_match_all(h, "src=\"data:")))
    for (raw in c(H, "<b>g1</b>", "<script>t", "<i>G01", "<script>x")) {
      expect_false(grepl(raw, h, fixed = TRUE), label = raw)
    }
    rel <- grep("^#", .rp_hrefs(h), value = TRUE, invert = TRUE)
    expect_true(all(file.exists(file.path(p, rel))))
    expect_false(any(grepl("^(/|[a-z]+:)", rel)))
    style <- .rp_match_all(h, "(?s)<style>.*?</style>")
    expect_false(any(grepl("@import|url\\(", style, ignore.case = TRUE)))
    ids <- .rp_ids(h)
    expect_false(anyDuplicated(ids) > 0)
  }
  hon <- .rp_html(on_img)
  expect_identical(length(.rp_match_all(hon, "src=\"data:image/png;base64,")),
                   if (capabilities("png")) 2L else 0L)
  expect_identical(length(.rp_match_all(.rp_html(off_img), "<img")), 0L)
  expect_true(grepl("&lt;script&gt;alert(1)&lt;/script&gt;", hon, fixed = TRUE))
  expect_true(grepl("a&amp;b &quot;c&quot; &#39;d&#39;", hon, fixed = TRUE))
  expect_true(grepl("&lt;b&gt;g1&lt;/b&gt;", hon, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;t&lt;/script&gt; &amp; co", hon, fixed = TRUE))
  hk <- .rp_html(pk)
  expect_true(grepl("&lt;i&gt;G01&lt;/i&gt;", hk, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;x", hk, fixed = TRUE))
})

test_that("inputs are checked and folders protected before anything is written (AC-007: folder input reporting differently, user files deleted, a table reported as provenance, themes of other data, partial folders)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  q1 <- q$q1; q2 <- q$q2
  x <- list(q1 = q1, q2 = q2)
  f1 <- file.path(base, "f1"); f2 <- file.path(base, "f2")
  writeRich(q1, f1); writeRich(q2, f2)
  pth <- function(n) file.path(base, n)

  p <- pth("p")
  richReport(x, p)
  richReport(c(q1 = f1, q2 = f2), pth("p3"))
  for (f in list.files(file.path(p, "tables"))) {
    expect_identical(unname(tools::md5sum(file.path(pth("p3"), "tables", f))),
                     unname(tools::md5sum(file.path(p, "tables", f))), label = f)
  }
  expect_setequal(list.files(file.path(pth("p3"), "tables")), list.files(file.path(p, "tables")))

  richReport(q1, pth("p4"))
  expect_identical(vapply(.rp_manifest(pth("p4"))$groups, `[[`, "", "name"), "Group 1")

  expect_error(richReport(bindRich(x), pth("p5")), "data frame.*analysis record")
  expect_error(richReport(c(f1, f1), pth("p6")), basename(f1), fixed = TRUE)

  d <- pth("d"); dir.create(d); writeLines("mine", file.path(d, "keep.txt"))
  before_d <- .rp_md5(d)
  expect_error(richReport(x, d), "not empty")
  expect_error(richReport(x, d, overwrite = TRUE), "not written by richReport")
  expect_identical(.rp_files(d), "keep.txt")
  expect_identical(.rp_md5(d), before_d)

  writeLines("my notes", file.path(p, "notes.txt"))
  before_p <- .rp_md5(p)
  expect_error(richReport(x, p), "not empty")
  expect_identical(.rp_md5(p), before_p)
  richReport(x, p, overwrite = TRUE)
  expect_identical(readLines(file.path(p, "notes.txt")), "my notes")
  expect_true(file.exists(file.path(p, "report.html")))

  qm <- .rp_fqm()
  expect_error(richReport(x, pth("p7"), themes = richThemes(list(q1 = q1, q2 = qm))), "themes")
  expect_error(richReport(x, pth("p8"), stats = .rp_st()), "stats")
  expect_error(richReport(list(q1 = q1, n = .rp_fnew(q1)), pth("p9")), "@detail", fixed = TRUE)
  expect_error(richReport(x, pth("p10"), figure.top = 0), "figure.top", fixed = TRUE)
  th <- richThemes(x)
  mem <- th$members$Annot[th$members$theme_id == "theme_001"]
  other <- setdiff(unique(mem), th$themes$representative[1])
  expect_gt(length(other), 0)
  th$themes$representative[1] <- other[1]
  expect_error(richReport(x, pth("p11"), themes = th), "themes")
  for (n in c("p5", "p6", "p7", "p8", "p9", "p10", "p11")) {
    expect_false(file.exists(pth(n)), label = n)
  }
})

test_that("rerun.R regenerates the report byte for byte from the saved results (AC-008: an unreproducible rerun, the failed group or labels lost, a rebuilt theme silently differing, a tampered table regenerated)", {
  skip_on_os("windows")
  base <- .rp_tmp()
  p <- file.path(base, "p")
  .rp_main(p)
  q <- file.path(base, "q")
  owd <- setwd(p)
  status <- system2(.rp_rscript(), c("rerun.R", shQuote(q)), env = .rp_lib_env(),
                    stdout = FALSE, stderr = FALSE)
  setwd(owd)
  expect_identical(status, 0L)
  .rp_expect_rerun_identical(p, q)
  expect_identical(.rp_decode(q, "themes")$label[1], "Edited theme")

  pc <- .rp_copy(p, file.path(base, "pc"))
  f <- file.path(pc, "tables/comparison.tsv")
  bytes <- readBin(f, "raw", file.info(f)$size)
  bytes[length(bytes) - 1L] <- as.raw(bitwXor(as.integer(bytes[length(bytes) - 1L]), 1L))
  writeBin(bytes, f)
  expect_error(richReport(pc, file.path(base, "q2")), "tables/comparison.tsv", fixed = TRUE)

  pt <- .rp_copy(p, file.path(base, "pt"))
  tf <- file.path(pt, "tables/themes.tsv")
  old_sha <- digest::digest(file = tf, algo = "sha256")
  lines <- readLines(tf, encoding = "UTF-8")
  head <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  rep_col <- match("representative", head)
  row <- which(startsWith(lines, "theme_001\t"))
  cells <- strsplit(lines[row], "\t", fixed = TRUE)[[1]]
  expect_identical(cells[rep_col], "K")
  cells[rep_col] <- "M"
  lines[row] <- paste(cells, collapse = "\t")
  con <- file(tf, "wb"); writeLines(lines, con, sep = "\n"); close(con)
  mf <- file.path(pt, "manifest.json")
  mtext <- readLines(mf, encoding = "UTF-8")
  expect_identical(sum(grepl(old_sha, mtext, fixed = TRUE)), 1L)
  mtext <- sub(old_sha, digest::digest(file = tf, algo = "sha256"), mtext, fixed = TRUE)
  con <- file(mf, "wb"); writeLines(mtext, con, sep = "\n"); close(con)
  expect_error(richReport(pt, file.path(base, "q3")), "tables/themes.tsv", fixed = TRUE)

  expect_error(richReport(p, file.path(base, "q4"), title = "x"), "title")
  for (n in c("q2", "q3", "q4")) expect_false(file.exists(file.path(base, n)), label = n)
})

test_that("a rerun under another session collation reproduces the figures (AC-008 q5: figures following the caller's LC_COLLATE)", {
  skip_on_os("windows")
  probe <- suppressWarnings(system2("env", c("-u", "LC_ALL", "LC_COLLATE=C.UTF-8", .rp_rscript(),
                                             "--vanilla", "-e",
                                             shQuote("cat(Sys.getlocale('LC_COLLATE'))")),
                                    stdout = TRUE, stderr = FALSE))
  skip_if_not(identical(probe, "C.UTF-8"), "the C.UTF-8 collation is unavailable")
  base <- .rp_tmp()
  p <- file.path(base, "p")
  .rp_main(p)
  q5 <- file.path(base, "q5")
  owd <- setwd(p)
  status <- system2("env", c("-u", "LC_ALL", "LC_COLLATE=C.UTF-8", .rp_lib_env(), .rp_rscript(),
                             "rerun.R", shQuote(q5)), stdout = FALSE, stderr = FALSE)
  setwd(owd)
  expect_identical(status, 0L)
  .rp_expect_rerun_identical(p, q5)
})

# ---- fix round 1: Codex M6.3 review and sync-audit test debt ----------------

test_that("overwrite removes the report's own entries literally and stops when one cannot be removed (Codex M6.3 #1: a bracketed folder name deleting a sibling folder; a failed removal leaving stale files)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  x <- list(q1 = q$q1, q2 = q$q2)
  p <- file.path(base, "study[1]")
  richReport(x, p, figures = FALSE)
  keep <- file.path(base, "study1", "results", "keep.txt")
  dir.create(dirname(keep), recursive = TRUE)
  writeLines("precious", keep)
  writeLines("stale", file.path(p, "results", "stale.txt"))
  richReport(x, p, figures = FALSE, overwrite = TRUE)
  expect_true(file.exists(keep))
  expect_false(file.exists(file.path(p, "results", "stale.txt")))
  expect_true(.rp_has(.rp_html(p), "Group 2: q2"))

  skip_on_os("windows")
  pl <- file.path(base, "pl")
  richReport(x, pl, figures = FALSE)
  locked <- file.path(pl, "tables", "locked")
  dir.create(locked)
  writeLines("x", file.path(locked, "a.txt"))
  writeLines("x", file.path(locked, "b.txt"))
  Sys.chmod(locked, "0555")
  on.exit(Sys.chmod(locked, "0755"), add = TRUE, after = FALSE)
  if (unlink(file.path(locked, "a.txt")) == 0L) skip("file removal is not refused here (root?)")
  expect_error(richReport(x, pl, figures = FALSE, overwrite = TRUE), "could not remove",
               fixed = TRUE)
  Sys.chmod(locked, "0755")
  richReport(x, pl, figures = FALSE, overwrite = TRUE)
  expect_false(dir.exists(locked))
})

test_that("figures do not depend on the session's graphics settings and leave them as found (Codex M6.3 #2: a rerun drawing other figures after update_geom_defaults() or another bitmapType)", {
  skip_on_os("windows")
  skip_if_not(capabilities("png"), "the png device is unavailable")
  skip_if_not(capabilities("cairo"), "cairo is unavailable")
  base <- .rp_tmp()
  old <- ggplot2::update_geom_defaults("point", list(shape = 8, fill = "red", stroke = 4))
  on.exit(ggplot2::update_geom_defaults("point", old), add = TRUE)
  bar <- ggplot2::update_geom_defaults("bar", list(colour = "red", linewidth = 3, alpha = 0.1))
  on.exit(ggplot2::update_geom_defaults("bar", bar), add = TRUE)
  op <- options(bitmapType = "cairo-png")
  on.exit(options(op), add = TRUE)
  p <- file.path(base, "p")
  .rp_main(p)
  expect_identical(ggplot2::GeomPoint$default_aes$shape, 8)
  expect_identical(getOption("bitmapType"), "cairo-png")
  q <- file.path(base, "q")
  owd <- setwd(p)
  status <- system2(.rp_rscript(), c("rerun.R", shQuote(q)), env = .rp_lib_env(),
                    stdout = FALSE, stderr = FALSE)
  setwd(owd)
  expect_identical(status, 0L)
  .rp_expect_rerun_identical(p, q)
  pg <- file.path(base, "pg"); qg <- file.path(base, "qg")
  richReport(.th_fgsea(FALSE), pg)
  owd <- setwd(pg)
  status <- system2(.rp_rscript(), c("rerun.R", shQuote(qg)), env = .rp_lib_env(),
                    stdout = FALSE, stderr = FALSE)
  setwd(owd)
  expect_identical(status, 0L)
  .rp_expect_rerun_identical(pg, qg)
  expect_identical(ggplot2::GeomBar$default_aes$alpha, 0.1)
})

test_that("themes whose settings were edited are refused (Codex M6.3 #3: a report showing settings its rerun does not reproduce)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  x <- list(q1 = q$q1, q2 = q$q2)
  edits <- list(background_size = 999999L, method = "jaccard", linkage = "complete")
  for (e in names(edits)) {
    th <- richThemes(x)
    expect_false(identical(th$settings[[e]], edits[[e]]), label = e)
    th$settings[[e]] <- edits[[e]]
    pe <- file.path(base, e)
    expect_error(richReport(x, pe, themes = th, figures = FALSE), "settings", label = e)
    expect_false(file.exists(pe), label = e)
  }
})

test_that("comparison cells come from their own group when names contain the key separator (Codex M6.3 #4: another group's state and adjusted p shown)", {
  base <- .rp_tmp()
  leg <- .rp_fleg(.rp_fq()$q1)
  two <- function(ids) {
    o <- leg
    o@result <- o@result[1:2, ]
    o@result$Annot <- ids
    o
  }
  x <- stats::setNames(list(two(c("B\001C", "zz")), two(c("yy", "C"))), c("A", "A\001B"))
  b <- bindRich(x)
  pick <- function(g, id) b[as.character(b$group) == g & b$Annot == id, ]
  expect_identical(c(pick("A", "C")$state, pick("A\001B", "C")$state),
                   c("unavailable", "reported"))
  expect_false(isTRUE(all.equal(pick("A\001B", "C")$Padj, pick("A", "B\001C")$Padj)))
  p <- file.path(base, "p")
  richReport(x, p, figures = FALSE)
  rows <- .rp_table_rows(.rp_section(.rp_html(p), "comparison"))
  cell <- Filter(function(r) identical(r[1], "C"), rows[-1])
  expect_length(cell, 1L)
  expect_identical(cell[[1]][-1],
                   c("unavailable", sprintf("reported (%.4g)", pick("A\001B", "C")$Padj)))
})

test_that("a result without an adjusted p column is compared with its adjusted p shown as NA (Codex M6.3 #5: comparison rendering crashing)", {
  base <- .rp_tmp()
  leg <- .rp_fleg(.rp_fq()$q1)
  leg@result$Padj <- NULL
  x <- list(a = leg, b = leg)
  expect_false("Padj" %in% names(bindRich(x)))
  p <- file.path(base, "p")
  richReport(x, p)
  rows <- .rp_table_rows(.rp_section(.rp_html(p), "comparison"))
  expect_identical(vapply(rows[-1], `[`, "", 1), sort(leg@result$Annot, method = "radix"))
  for (r in rows[-1]) expect_identical(r[-1], rep("reported (NA)", 2), label = r[1])
  if (isTRUE(unname(capabilities("png")))) {
    expect_length(.rp_match_all(.rp_html(p), "<img "), 2L)
  }
  q <- file.path(base, "q")
  richReport(p, q)
  expect_identical(.rp_files(p), .rp_files(q))
  fixed <- grep("^(tables/|results/|rerun\\.R$)", .rp_files(p), value = TRUE)
  expect_gt(length(fixed), 0)
  expect_identical(unname(tools::md5sum(file.path(p, fixed))),
                   unname(tools::md5sum(file.path(q, fixed))))
  expect_identical(.rp_match_all(.rp_html(p), "<img [^>]*>"),
                   .rp_match_all(.rp_html(q), "<img [^>]*>"))
})

test_that("numbers are shown in full when whole, with 4 significant digits otherwise, and NA when missing (REQ-011; F1: a large count shown as 1.235e+05)", {
  expect_identical(richR:::.report_num(c(123456L, NA)), c("123456", "NA"))
  expect_identical(richR:::.report_num(c(123456, 1e15, -2e5, 0, 0.1, 2 / 3, NA)),
                   c("123456", "1000000000000000", "-200000", "0", "0.1", "0.6667", "NA"))
})

test_that("the generation time is UTC whatever the session time zone (REQ-005; F3: local time stamped Z)", {
  base <- .rp_tmp()
  old <- Sys.getenv("TZ", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("TZ") else Sys.setenv(TZ = old), add = TRUE)
  Sys.setenv(TZ = "Asia/Tokyo")
  p <- file.path(base, "p")
  richReport(.rp_fq()$q1, p, figures = FALSE)
  generated <- .rp_manifest(p)$generated
  at <- as.POSIXct(generated, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  expect_lt(abs(as.numeric(difftime(at, Sys.time(), units = "mins"))), 5)
  expect_true(.rp_has(.rp_html(p), paste("Generated:", generated)))
})

test_that("identifier text and record values are escaped at every site (AC-006; F5: a hostile identifier raw in comparison rows, member links or the representative, a hostile id raw in a record value)", {
  base <- .rp_tmp()
  hid <- "<K&'\">"
  gene <- "<i>G01</i>"
  a <- .rp_a80(c(G01 = gene))
  a$Term[a$Term == "K"] <- hid
  query <- replace(.rp_G(c(1:6, 11:14)), 1, gene)
  x <- list(a = .rp_enrich(query, annot = a), b = .rp_enrich(.rp_G(c(1:8, 30:35)), annot = a))
  expect_true(hid %in% x$a@result$Annot)
  th <- richThemes(x)
  expect_identical(th$themes$representative[1], hid)
  p <- file.path(base, "p")
  richReport(x, p, themes = th, figures = FALSE)
  h <- .rp_html(p)
  for (raw in c("<K&", "<i>G01")) expect_false(grepl(raw, h, fixed = TRUE), label = raw)
  rows <- .rp_table_rows(.rp_section(h, "comparison"))
  expect_true(hid %in% vapply(rows[-1], `[`, "", 1))
  t1 <- .rp_section(h, "theme_001", until = "<h3")
  expect_true(.rp_has(t1, paste0("Representative: ", hid)))
  expect_true(grepl(">&lt;K&amp;&#39;&quot;&gt;</a>", t1, fixed = TRUE))
  expect_identical(.rp_record_value(.rp_section(h, "group-1"), "input$ids"),
                   paste(query, collapse = ", "))
})

test_that("a rerun reads each group only from its own results folder (F7: an edited manifest.json reading results from elsewhere)", {
  base <- .rp_tmp()
  q <- .rp_fq()
  p <- file.path(base, "p")
  richReport(list(q1 = q$q1, q2 = q$q2), p, figures = FALSE)
  writeRich(q$q2, file.path(base, "ext"))
  mtext <- readLines(file.path(p, "manifest.json"), encoding = "UTF-8")
  line <- "\"folder\": \"results/group_1\""
  expect_identical(sum(grepl(line, mtext, fixed = TRUE)), 1L)
  for (k in 1:2) {
    to <- c("../ext", "results/group_2")[k]
    pc <- .rp_copy(p, file.path(base, paste0("pc", k)))
    con <- file(file.path(pc, "manifest.json"), "wb")
    writeLines(sub(line, sprintf("\"folder\": \"%s\"", to), mtext, fixed = TRUE), con, sep = "\n")
    close(con)
    out <- file.path(base, paste0("q", k))
    expect_error(richReport(pc, out), "results/group_1", fixed = TRUE, label = to)
    expect_false(file.exists(out), label = to)
  }
})

# Integration seam: SPEC-RICHR-009's new status must survive SPEC-RICHR-008.
test_that("reports and reruns retain undetermined zero-overlap theme evidence (SPEC-RICHR-008/009)", {
  base <- .rp_tmp()
  annot <- data.frame(GeneID = paste0("G", 1:60),
                       Annot = rep(paste0("T", 1:6), each = 10))
  annot$Term <- annot$Annot
  run <- function(q, ...) suppressMessages(enrich(paste0("G", q), annot,
    pvalue = 1, padj = 1, minSize = 1, minGSSize = 1, ...))
  groups <- list(reported = run(c(1:8, 11:12)),
                 all_family = run(c(11:18, 21:24)),
                 overlap_family = run(c(11:18, 21:24)))
  # A supported historical record can include zero-overlap tests even though
  # the current upstream engine does not create that policy.
  m <- enrichmentInfo(groups$all_family)
  m$testing$zero_overlap_count <- 1L
  m$testing$tested_count <- as.integer(nrow(groups$all_family@tested) + 1L)
  m$testing$family_policy <- "explicit_historical_fixture"
  attr(groups$all_family, "richR_metadata") <- m
  bound <- bindRich(groups)
  expect_identical(bound$state[bound$Annot == "T1"],
                   c("reported", "untested_or_zero_overlap", "not_tested"))
  th <- richThemes(groups)
  tid <- th$members$theme_id[th$members$Annot == "T1"][1]
  view <- themeGenes(th)
  rows <- view[view$theme_id == tid, ]
  expect_identical(unique(rows$group_status[rows$group == "all_family"]), "undetermined")
  expect_identical(unique(rows$group_status[rows$group == "overlap_family"]), "absent")

  p <- file.path(base, "p"); q <- file.path(base, "q")
  richReport(groups, p, themes = th, figures = FALSE)
  richReport(p, q)
  for (folder in c(p, q)) {
    expect_identical(.rp_decode(folder, "theme_gene_view"), view)
    html <- .rp_section(.rp_html(folder), tid, until = "<h3")
    table <- .rp_table_rows(html)
    cols <- match(c("group", "group_status"), table[[1]])
    all_rows <- Filter(function(r) identical(r[cols[1]], "all_family"), table[-1])
    expect_length(all_rows, sum(rows$group == "all_family"))
    expect_gt(length(all_rows), 0)
    expect_identical(unique(vapply(all_rows, `[`, "", cols[2])), "undetermined")
  }
  .rp_expect_rerun_identical(p, q)
})

# Literal deletion must preserve R's usual home-directory expansion.
test_that("report overwrite expands a tilde path once and returns the supplied path", {
  home <- path.expand("~")
  base <- tempfile("richReport-tilde-", tmpdir = home)
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE, expand = FALSE), add = TRUE)
  path <- file.path("~", basename(base), "study[1]")
  disk <- path.expand(path)
  x <- .rp_fq()$q1
  expect_identical(richReport(x, path, figures = FALSE), path)
  writeLines("stale", file.path(disk, "results", "stale.txt"))
  keep <- file.path(base, "study1", "results", "keep.txt")
  dir.create(dirname(keep), recursive = TRUE)
  writeLines("precious", keep)
  expect_identical(richReport(x, path, figures = FALSE, overwrite = TRUE), path)
  expect_false(file.exists(file.path(disk, "results", "stale.txt")))
  expect_identical(readLines(keep), "precious")
  expect_identical(readRich(file.path(disk, "results/group_1")), x)
})

# Accepted legacy tables need not contain every statistic a figure requires.
for (removed in list("Significant", "Annotated", c("Significant", "Annotated"))) local({
  columns <- removed
  test_that(paste("ORA report keeps legacy data without plotting columns", paste(columns, collapse = ", ")), {
    base <- .rp_tmp()
    obj <- .rp_fleg(.rp_fq()$q1)
    obj@result[columns] <- NULL
    input <- file.path(base, "input")
    writeRich(obj, input)
    expect_identical(readRich(input), obj)
    p <- file.path(base, "p"); q <- file.path(base, "q")
    richReport(obj, p)
    expect_identical(readRich(file.path(p, "results/group_1")), obj)
    expect_identical(.rp_decode(p, "comparison"), bindRich(list(obj)))
    expected <- if (!isTRUE(unname(capabilities("png")))) {
      "No figure: the png device is unavailable."
    } else paste0("No figure: required plotting columns missing: ", paste(columns, collapse = ", "), ".")
    expect_identical(.rp_captions(.rp_html(p)), expected)
    expect_length(.rp_match_all(.rp_html(p), "<img "), 0L)
    richReport(p, q)
    expect_identical(.rp_captions(.rp_html(q)), expected)
    expect_identical(readRich(file.path(q, "results/group_1")), obj)
    expect_identical(.rp_decode(q, "comparison"), bindRich(list(obj)))
  })
})

test_that("GSEA report keeps accepted results without an NES column and omits the figure", {
  base <- .rp_tmp()
  obj <- .th_fgsea(FALSE)$a
  obj@result$NES <- NULL
  input <- file.path(base, "input")
  writeRich(obj, input)
  expect_identical(readRich(input), obj)
  p <- file.path(base, "p"); q <- file.path(base, "q")
  richReport(obj, p)
  expect_identical(readRich(file.path(p, "results/group_1")), obj)
  expect_identical(.rp_decode(p, "comparison"), bindRich(list(obj)))
  expected <- if (!isTRUE(unname(capabilities("png")))) {
    "No figure: the png device is unavailable."
  } else "No figure: required plotting columns missing: NES."
  expect_identical(.rp_captions(.rp_html(p)), expected)
  expect_length(.rp_match_all(.rp_html(p), "<img "), 0L)
  richReport(p, q)
  expect_identical(.rp_captions(.rp_html(q)), expected)
  expect_identical(readRich(file.path(q, "results/group_1")), obj)
  expect_identical(.rp_decode(q, "comparison"), bindRich(list(obj)))
})
