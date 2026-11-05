# compareGSEA() significance cutoffs.
#
# compareGSEA() keeps a group's pathways by the caller's cutoff, compared
# strictly (<): with a non-NULL `padj` the adjusted cutoff is used in place of
# the raw one; with `padj = NULL` the raw cutoff applies together with the
# `padj < 1` guard (the same rule compareResult() documents for ORA results).
# The cutoffs must never be resolved to the same-named columns of the result
# table (`padj`, `pval`).
#
# Fixtures are built inline; expectations are valid in testthat editions 2 and 3.

cutoffs_scores <- function() {
  genes <- sprintf("g%02d", 1:24)
  stats::setNames(seq(2.3, -2.3, length.out = length(genes)), genes)
}

cutoffs_annot <- function() {
  genes <- names(cutoffs_scores())
  data.frame(
    GeneID = genes,
    Term = rep(c("P1", "P2", "P3", "P4", "P5", "P6"), each = 4),
    stringsAsFactors = FALSE
  )
}

# One group with six pathways whose raw and adjusted p-values separate the
# cutoff cases:
#   P1 passes both cutoffs; P2 passes the raw cutoff only; P3 passes the
#   adjusted cutoff only; P4 passes both; P5 passes neither; P6 has a small raw
#   p-value and an adjusted p-value of exactly 1.
cutoffs_result <- function() {
  fc <- cutoffs_scores()
  res <- data.frame(
    pathway = c("P1", "P2", "P3", "P4", "P5", "P6"),
    pval = c(0.01, 0.03, 0.06, 0.04, 0.20, 0.02),
    padj = c(0.02, 0.09, 0.04, 0.07, 0.30, 1.00),
    log2err = rep(0.5, 6),
    ES = c(0.9, 0.8, -0.7, 0.6, 0.3, -0.5),
    NES = c(1.9, 1.7, -1.5, 1.3, 0.6, -1.1),
    size = rep(4L, 6),
    leadingEdge = rep("g01", 6),
    stringsAsFactors = FALSE
  )
  methods::new(
    "GSEAResult", result = res, pvalueCutoff = 1, pAdjustMethod = "BH",
    padjCutoff = numeric(), genenumber = length(fc), organism = character(),
    gene = names(fc), input = fc, ontology = character(), keytype = "",
    sep = ","
  )
}

# The set of pathways drawn for each group in a compareGSEA() plot.
cutoffs_drawn <- function(plot_obj) {
  sort(unique(as.character(plot_obj$data$Group)))
}

test_that("compareGSEA() with a padj cutoff keeps pathways passing the adjusted cutoff", {
  # Catches: a padj cutoff stopping with "object 'Padj' not found".
  skip_if_not_installed("fgsea")
  p <- compareGSEA(list(A = cutoffs_result()), cutoffs_annot(),
                   pvalue = 0.05, padj = 0.08)
  expect_s3_class(p, "ggplot")
  expect_identical(cutoffs_drawn(p), c("P1", "P3", "P4"))
})

test_that("compareGSEA() with padj = NULL applies the raw cutoff and the padj < 1 guard", {
  skip_if_not_installed("fgsea")
  p <- compareGSEA(list(A = cutoffs_result()), cutoffs_annot(),
                   pvalue = 0.05, padj = NULL)
  expect_s3_class(p, "ggplot")
  expect_identical(cutoffs_drawn(p), c("P1", "P2", "P4"))
})

test_that("compareGSEA() compares cutoffs strictly at the boundary", {
  skip_if_not_installed("fgsea")
  # P4 has pval exactly 0.04; P3 has padj exactly 0.04.
  by_raw <- compareGSEA(list(A = cutoffs_result()), cutoffs_annot(),
                        pvalue = 0.04, padj = NULL)
  expect_identical(cutoffs_drawn(by_raw), c("P1", "P2"))
  by_adj <- compareGSEA(list(A = cutoffs_result()), cutoffs_annot(),
                        pvalue = 0.1, padj = 0.04)
  expect_identical(cutoffs_drawn(by_adj), "P1")
})

test_that("compareGSEA() uses the caller's padj, not the result's padj column", {
  # Catches: the cutoff resolving to the result's own padj column.
  skip_if_not_installed("fgsea")
  # With padj = 1 every pathway except P6 (padj exactly 1) passes.
  p <- compareGSEA(list(A = cutoffs_result()), cutoffs_annot(),
                   pvalue = 0.05, padj = 1)
  expect_identical(cutoffs_drawn(p), c("P1", "P2", "P3", "P4", "P5"))
})
