# ORA plot display labels. The results are built by hand, only exported plot
# functions are called, and the expectations are valid in testthat editions 2
# and 3.

plot_result <- function(ids, terms, level2 = "L2") {
  n <- length(ids)
  methods::new(
    "richResult",
    result = data.frame(
      Annot = ids, Term = terms,
      Annotated = rep(20, n), Significant = rep(5, n) + seq_len(n),
      RichFactor = (rep(5, n) + seq_len(n)) / 20,
      FoldEnrichment = rep(2, n) + seq_len(n) / 10,
      zscore = rep(3, n), Pvalue = 0.001 * seq_len(n), Padj = 0.002 * seq_len(n),
      GeneID = vapply(seq_len(n), function(i) {
        paste(sprintf("g%02d", seq(i, i + 5 + i)), collapse = ",")
      }, character(1)),
      Level2 = rep(level2, n),
      stringsAsFactors = FALSE, row.names = ids
    ),
    detail = data.frame(TERM = character(), Annot = character(), GeneID = character(),
                        Pvalue = numeric(), Padj = numeric(), stringsAsFactors = FALSE),
    pvalueCutoff = 0.05, pAdjustMethod = "BH", padjCutoff = 0.05, genenumber = 30,
    organism = "custom", ontology = "", gene = sprintf("g%02d", 1:30), keytype = "",
    sep = ","
  )
}

# Every category string a built ggplot shows: discrete axis labels, polar
# angle labels and text-layer labels.
plot_categories <- function(p) {
  b <- ggplot2::ggplot_build(p)
  pp <- b$layout$panel_params[[1]]
  out <- character()
  for (ax in c("x", "y")) {
    sc <- pp[[ax]]
    if (!is.null(sc) && is.function(sc$get_labels)) {
      out <- c(out, tryCatch(as.character(sc$get_labels()), error = function(e) character()))
    }
  }
  if (!is.null(pp$theta.labels)) out <- c(out, as.character(pp$theta.labels))
  for (d in b$data) if ("label" %in% names(d)) out <- c(out, as.character(d$label))
  unique(out[!is.na(out)])
}

test_that("richDot() shortens long labels when short = TRUE", {
  # Catches: richDot() accepting short but never passing it on.
  res <- plot_result(c("W1", "W2"), c("a b c d e f g h", "other"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (object in list(res, res@result)) {
    long <- grep("^a", plot_categories(richDot(object, pvalue = 1, padj = 1, short = TRUE)),
                 value = TRUE)
    expect_equal(length(long), 1)
    expect_true(grepl("\n", long[1], fixed = TRUE))
    plain <- grep("^a", plot_categories(richDot(object, pvalue = 1, padj = 1)), value = TRUE)
    expect_identical(plain, "a b c d e f g h")
  }
})
