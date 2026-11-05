# compareGSEA() per-group value tables.
#
# compareGSEA() must key its NES / adjusted-p / raw-p rows by group and pathway
# identifier, so groups whose retained pathway sets or orders differ still get
# each value from the matching (group, pathway) row of their own result.
#
# Fixtures are built inline; expectations are valid in testthat editions 2 and 3.

rows_scores <- function() {
  genes <- sprintf("g%02d", 1:20)
  stats::setNames(seq(2.1, -2.1, length.out = length(genes)), genes)
}

rows_annot <- function() {
  data.frame(
    GeneID = names(rows_scores()),
    Term = rep(c("P1", "P2", "P3", "P4", "P5"), each = 4),
    stringsAsFactors = FALSE
  )
}

rows_gsea <- function(pathway, pval, NES) {
  fc <- rows_scores()
  n <- length(pathway)
  res <- data.frame(
    pathway = pathway, pval = pval, padj = pmin(1, pval * 2),
    log2err = rep(0.5, n), ES = NES / 2, NES = NES, size = rep(4L, n),
    leadingEdge = rep("g01", n), stringsAsFactors = FALSE
  )
  methods::new(
    "GSEAResult", result = res, pvalueCutoff = 1, pAdjustMethod = "BH",
    padjCutoff = numeric(), genenumber = length(fc), organism = character(),
    gene = names(fc), input = fc, ontology = character(), keytype = "",
    sep = ","
  )
}

# Group A retains P4, P1, P2 (in that row order); group B retains P3, P2.
# P5 fails the raw cutoff in both groups.
rows_groups <- function() {
  list(
    A = rows_gsea(c("P4", "P1", "P2", "P5"), c(0.001, 0.002, 0.003, 0.5),
                  c(1.4, 1.1, 1.2, 0.3)),
    B = rows_gsea(c("P3", "P2", "P5", "P1"), c(0.001, 0.004, 0.6, 0.7),
                  c(-2.3, 2.2, 0.2, -0.1))
  )
}

rows_label_data <- function(plot_obj) {
  is_text <- vapply(plot_obj$layers, function(l) inherits(l$geom, "GeomText"),
                    logical(1))
  as.data.frame(plot_obj$layers[[which(is_text)[1]]]$data)
}

test_that("compareGSEA() labels every drawn (group, pathway) pair exactly once", {
  skip_if_not_installed("fgsea")
  p <- compareGSEA(rows_groups(), rows_annot(), pvalue = 0.05, label = TRUE)
  expect_s3_class(p, "ggplot")
  drawn <- sort(unique(paste(p$data$group, p$data$Group, sep = "@")))
  expect_identical(drawn, c("A@P1", "A@P2", "A@P4", "B@P2", "B@P3"))
  labels <- rows_label_data(p)
  labelled <- sort(paste(labels$group, labels$Group, sep = "@"))
  expect_identical(labelled, drawn)
})

test_that("compareGSEA() takes each NES from the matching group and pathway row", {
  skip_if_not_installed("fgsea")
  groups <- rows_groups()
  p <- compareGSEA(groups, rows_annot(), pvalue = 0.05, label = TRUE)
  labels <- rows_label_data(p)
  expected <- vapply(seq_len(nrow(labels)), function(i) {
    res <- groups[[labels$group[i]]]@result
    res$NES[res$pathway == labels$Group[i]]
  }, numeric(1))
  expect_equal(labels$NES, expected)
  expect_identical(labels$annotateText,
                   paste0(labels$group, ":(NES=", round(expected, 2), ")\n"))
})

test_that("compareGSEA() keeps rows aligned when one group retains a single pathway", {
  skip_if_not_installed("fgsea")
  groups <- rows_groups()
  groups$B <- rows_gsea(c("P5", "P3"), c(0.9, 0.01), c(0.1, -1.7))
  p <- compareGSEA(groups, rows_annot(), pvalue = 0.05, label = TRUE)
  labels <- rows_label_data(p)
  b_rows <- labels[labels$group == "B", , drop = FALSE]
  expect_identical(as.character(b_rows$Group), "P3")
  expect_equal(b_rows$NES, -1.7)
})
