# Shared fixtures of the theme tests (SPEC-RICHR-007), built from base R,
# testthat and the public API, with no network.
#
# F-TH   two hand-built ORA groups (no tested table): g1 reports A..F, g2
#        reports D (another Term) and G (Term NA). Kappa A-B = A-C = D-E =
#        0.8224, B-C = 0.6447, every other pair 0, on a 27-gene background.
# F-LINK one hand-built group whose average-linkage partition at 0.5 differs
#        from both single and complete linkage.
# F-TH2  two enrich() runs over one annotation: K and M share their overlap
#        in q1, K is tested-not-reported in q2, W is tested only.
# F-ZERO an enrich() run that reports nothing but tests Z.
# F-LEG  an as(<data.frame>, "richResult") object without a tested table.
# F-GSEA the characterization GSEA fixtures, in both leadingEdge forms.
#
# Every builder is paired with a role check; tests assert the roles before
# relying on them, so a fixture drift fails loudly instead of silently
# changing what a test checks.

.th_g <- function(i) sprintf("g%02d", i)

.th_res <- function(ids, sets, padj, terms) {
  df <- data.frame(Annot = ids, Term = terms, Annotated = lengths(sets) + 2L,
                   Significant = lengths(sets), Pvalue = padj, Padj = padj,
                   GeneID = vapply(sets, paste, "", collapse = ","),
                   stringsAsFactors = FALSE)
  rownames(df) <- ids
  methods::new("richResult", result = df, sep = ",")
}

# The F-TH evidence sets, by identifier.
.th_fth_sets <- function() {
  g <- .th_g
  list(A = g(1:8), B = g(c(1:7, 9)), C = g(c(2:8, 10)), D = g(20:27),
       E = g(c(20:26, 28)), F = g(35:38), G = g(40:43))
}

.th_fth <- function() {
  s <- .th_fth_sets()
  g1 <- .th_res(c("A", "B", "C", "D", "E", "F"), s[c("A", "B", "C", "D", "E", "F")],
                c(0.04, 0.001, 0.01, 0.02, 0.03, 0.005),
                c("alpha", "beta", "gamma", "death receptor", "epsilon", "phi"))
  g2 <- .th_res(c("D", "G"), list(.th_g(20:24), s$G), c(0.002, 0.003),
                c("death receptor (g2 label)", NA))
  list(g1 = g1, g2 = g2)
}

# Cohen's kappa of two gene sets against a background, from the 2 x 2 table,
# computed here and never read from package code.
.th_kappa <- function(a, b, U) {
  n11 <- length(intersect(a, b)); n10 <- length(setdiff(a, b))
  n01 <- length(setdiff(b, a)); n00 <- length(setdiff(U, union(a, b)))
  N <- n11 + n10 + n01 + n00
  if (n11 == 0) return(0)
  po <- (n11 + n00) / N
  pe <- ((n11 + n10) * (n11 + n01) + (n01 + n00) * (n10 + n00)) / N^2
  if (pe == 1) 0 else (po - pe) / (1 - pe)
}

# Theme member sets as sorted strings ("ABC"), sorted, for order-free
# comparison of partitions.
.th_partition <- function(th, id = "Annot") {
  m <- unique(th$members[, c("theme_id", id)])
  unname(sort(vapply(split(m[[id]], m$theme_id),
                     function(v) paste(sort(v, method = "radix"), collapse = ""), ""),
              method = "radix"))
}

.th_fth_roles <- function() {
  s <- .th_fth_sets(); U <- unique(unlist(s))
  testthat::expect_identical(length(U), 27L)
  k <- function(a, b) .th_kappa(s[[a]], s[[b]], U)
  testthat::expect_equal(round(c(k("A", "B"), k("A", "C"), k("D", "E"), k("B", "C")), 4),
                         c(0.8224, 0.8224, 0.8224, 0.6447))
  ids <- names(s)
  for (a in ids) for (b in ids) {
    if (a < b && !paste0(a, b) %in% c("AB", "AC", "BC", "DE")) {
      testthat::expect_identical(k(a, b), 0)
    }
  }
}

.th_flink_sets <- function() {
  h <- function(i) sprintf("h%02d", i)
  list(P = h(1:10), Q = h(3:12), R = h(6:15), T = h(1:4), Z = h(40:79))
}

.th_flink <- function() {
  s <- .th_flink_sets()
  .th_res(names(s), s, c(0.01, 0.02, 0.03, 0.04, 0.05), paste("term", names(s)))
}

# F-LINK roles: at cutoff 0.5 average linkage joins R to {P,Q} (mean kappa
# 0.511 >= 0.5) and leaves T alone (mean 0.362 to {P,Q,R}); single linkage
# would also join T (P-T 0.5217), complete linkage would not join R (P-R
# 0.3889).
.th_flink_roles <- function() {
  s <- .th_flink_sets(); U <- unique(unlist(s))
  testthat::expect_identical(length(U), 55L)
  k <- function(a, b) .th_kappa(s[[a]], s[[b]], U)
  testthat::expect_equal(round(c(k("P", "Q"), k("Q", "R"), k("P", "T"), k("P", "R"),
                                 k("Q", "T"), k("R", "T")), 4),
                         c(0.7556, 0.6333, 0.5217, 0.3889, 0.2029, 0))
  testthat::expect_gte((k("P", "R") + k("Q", "R")) / 2, 0.5)
  testthat::expect_lt((k("P", "T") + k("Q", "T") + k("R", "T")) / 3, 0.5)
  testthat::expect_gte(k("P", "T"), 0.5)
  testthat::expect_lt(k("P", "R"), 0.5)
}

.th_fth2_annot <- function() {
  G <- sprintf("G%02d", 1:80)
  mk <- function(id, label, idx) {
    data.frame(GeneID = G[idx], Term = id, Annot = label, stringsAsFactors = FALSE)
  }
  rbind(mk("K", "kinase signalling", 1:8), mk("M", "kinase signalling (sub)", c(1:6, 9)),
        mk("L", "lipid handling", 11:16), mk("N", "nucleotide pathway", 30:37),
        mk("W", "weakly hit set", 17:29), mk("Z", "large background set", 41:80))
}

.th_enrich <- function(q, pvalue = 0.05, sep = ",") {
  suppressMessages(richR::enrich(q, .th_fth2_annot(), pvalue = pvalue, padj = NULL,
                                 minSize = 2, minGSSize = 3, sep = sep))
}

.th_fth2 <- function() {
  G <- sprintf("G%02d", 1:80)
  list(q1 = .th_enrich(G[c(1:6, 11:14, 17)]), q2 = .th_enrich(G[c(7, 8, 30:35)]))
}

.th_state <- function(b, id, group) {
  b$state[b$Annot == id & as.character(b$group) == group]
}

.th_fth2_roles <- function(q1, q2) {
  b <- richR::bindRich(list(q1 = q1, q2 = q2))
  st <- function(id, g) .th_state(b, id, g)
  testthat::expect_identical(
    c(st("K", "q1"), st("M", "q1"), st("L", "q1"), st("W", "q1"), st("N", "q1")),
    c("reported", "reported", "reported", "tested_not_reported", "not_tested"))
  testthat::expect_identical(
    c(st("N", "q2"), st("K", "q2"), st("M", "q2"), st("L", "q2"), st("W", "q2")),
    c("reported", "tested_not_reported", "not_tested", "not_tested", "not_tested"))
  testthat::expect_identical(b$GeneID[b$Annot == "K" & b$group == "q2"], "G07,G08")
  testthat::expect_identical(b$GeneID[b$Annot == "K" & b$group == "q1"],
                             "G01,G02,G03,G04,G05,G06")
  testthat::expect_false(identical(richR::enrichmentInfo(q1), richR::enrichmentInfo(q2)))
}

.th_fzero <- function() .th_enrich(c("G41", "G42"), pvalue = 1e-20)

.th_fleg <- function() {
  methods::as(data.frame(Annot = c("K", "P"), Term = c("kinase (legacy)", "P term"),
                         Annotated = c(8L, 2L), Significant = c(3L, 2L),
                         GeneID = c("G01,G02,G03", "G50,G51"), Pvalue = c(0.01, 0.02),
                         Padj = c(0.01, 0.02), stringsAsFactors = FALSE), "richResult")
}

# F-GSEA: both rankings in one leadingEdge form (TRUE = joined text, FALSE =
# list); every pathway reported.
.th_fgsea <- function(table) {
  run <- function(name) {
    .fgsea_run(name, pvalue = 1, padj = NULL, minSize = 3, maxSize = 20, table = table)
  }
  list(a = run("gsea_ranked_a.csv"), b = run("gsea_ranked_b.csv"))
}

# T5 has an empty leading edge in a GSEA result, in either storage form.
.th_t5_empty <- function(x) {
  le <- x@result$leadingEdge[x@result$pathway == "T5"]
  if (is.list(le)) length(le[[1]]) == 0L else identical(le, "")
}
