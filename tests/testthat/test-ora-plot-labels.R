# ORA plot display labels. The results are built by hand, only exported plot
# functions are called, and the expectations are valid in testthat editions 2
# and 3.
#
# Two identifiers sharing a label must be drawn as two categories
# ("label [ID]"), never merged and never an error.

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

shared_result <- function() plot_result(c("A1", "A2", "B1"), c("Shared", "Shared", "Unique"))

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

# Every text label drawn anywhere in a plot's grob tree (for composed plots,
# such as richClusterDot, whose panels are embedded grobs).
grob_labels <- function(p) {
  walk <- function(g) {
    out <- character()
    if (inherits(g, "text") && !is.null(g$label)) out <- as.character(g$label)
    kids <- c(if (inherits(g, "gTree")) g$children, if (!is.null(g$grobs)) g$grobs)
    for (k in kids) out <- c(out, walk(k))
    out
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  unique(walk(ggplot2::ggplotGrob(p)))
}

fold <- function(x) gsub("\\s+", " ", x)

# Evaluate one plot call; an error is returned (and reported) instead of
# aborting the remaining plot paths.
attempt <- function(expr) tryCatch(expr, error = function(e) e)

expect_one_category_per_id <- function(p, what) {
  expect_false(inherits(p, "error"),
               info = paste(what, if (inherits(p, "error")) conditionMessage(p)))
  if (inherits(p, "error")) return(invisible())
  expect_true(inherits(p, "ggplot"), info = what)
  cats <- fold(plot_categories(p))
  expect_true(all(c("Shared [A1]", "Shared [A2]", "Unique") %in% cats), info = what)
  expect_false("Shared" %in% cats, info = what)
}

fc <- stats::setNames(seq(2, -2, length.out = 30), sprintf("g%02d", 1:30))

test_that("every ORA plot draws one category per identifier when labels are shared", {
  # Catches: ORA plots that merge identifiers sharing a label, or error on them.
  res <- shared_result()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- list(
    richBar = attempt(richBar(res, pvalue = 1, padj = 1)),
    richBar_order = attempt(richBar(res, pvalue = 1, padj = 1, order = TRUE)),
    richDot = attempt(richDot(res, pvalue = 1, padj = 1)),
    richDot_order = attempt(richDot(res, pvalue = 1, padj = 1, order = TRUE)),
    richNetplot = attempt(richNetplot(res, pvalue = 1, padj = 1)),
    richLollipop = attempt(richLollipop(res, pvalue = 1, padj = 1)),
    richVolcano = attempt(richVolcano(res, pvalue = 1, padj = 1)),
    richCircle = attempt(richCircle(res, pvalue = 1, padj = 1)),
    richScatter = attempt(richScatter(res, pvalue = 1, padj = 1)),
    richTermSim = attempt(richTermSim(res, pvalue = 1, padj = 1, sim.cutoff = 0)),
    richGeneDot = attempt(richGeneDot(res, fc = fc, pvalue = 1, padj = 1)),
    richGeneHeat = attempt(richGeneHeat(res, fc = fc, pvalue = 1, padj = 1)),
    richGeneBar = attempt(richGeneBar(res, pvalue = 1, padj = 1))
  )
  for (what in names(plots)) expect_one_category_per_id(plots[[what]], what)
})

test_that("richNetwork names one node per identifier when labels are shared", {
  # Catches: duplicate row names when two identifiers share a label.
  res <- shared_result()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  # igraph >= 2.0 deprecation notices from richNetwork's igraph calls are not
  # what this test is about; only those are muffled.
  p <- withCallingHandlers(
    richNetwork(res, pvalue = 1, padj = 1, weightcut = 0),
    lifecycle_warning_deprecated = function(w) invokeRestart("muffleWarning")
  )
  expect_one_category_per_id(p, "richNetwork")
})

two_results <- function() {
  list(
    first = plot_result(c("A1", "B1"), c("Shared", "Unique")),
    second = plot_result(c("A2", "A1"), c("Shared", "Shared"))
  )
}

test_that("richHeatmap matches rows across results by identifier", {
  # Catches: a heatmap joining two results by label.
  rr <- two_results()
  data <- richHeatmap(rr, returnData = TRUE)
  expect_identical(sort(data$Term), c("Shared [A1]", "Shared [A2]", "Unique"))
  a1 <- data[data$Term == "Shared [A1]", 2:3]
  expect_equal(unname(unlist(a1)), c(rr$first@result["A1", "Pvalue"],
                                     rr$second@result["A1", "Pvalue"]))
  a2 <- data[data$Term == "Shared [A2]", 2:3]
  expect_equal(unname(unlist(a2)), c(1, rr$second@result["A2", "Pvalue"]))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_true(inherits(richHeatmap(rr), "ggplot"))
})

test_that("richCompareDot and richClusterDot keep one row per identifier across results", {
  # Catches: multi-result dot plots merging identifiers that share a label.
  compared <- compareResult(two_results())
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- richCompareDot(compared, pvalue = 1, padj = 1)
  expect_one_category_per_id(p, "richCompareDot")
  b <- ggplot2::ggplot_build(p)
  y_labels <- fold(b$layout$panel_params[[1]]$y$get_labels())
  points_y <- y_labels[b$data[[1]]$y]
  # A1 is in both results: one row holding two points; A2 only in the second.
  expect_equal(sum(points_y == "Shared [A1]"), 2)
  expect_equal(sum(points_y == "Shared [A2]"), 1)
  cl <- attempt(richClusterDot(compared))
  expect_false(inherits(cl, "error"), info = "richClusterDot")
  if (!inherits(cl, "error")) {
    cats <- fold(grob_labels(cl))
    expect_true(all(c("Shared [A1]", "Shared [A2]", "Unique") %in% cats))
    expect_false("Shared" %in% cats)
  }
})

test_that("short = TRUE keeps labels that differ only in whitespace apart", {
  # Catches: shortening re-merging two different labels into one category.
  res <- plot_result(c("W1", "W2", "W3"), c("a b c d e f g h", "a  b c d e f g h", "other"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- list(
    richBar = attempt(richBar(res, pvalue = 1, padj = 1, short = TRUE)),
    richDot = attempt(richDot(res, pvalue = 1, padj = 1, short = TRUE)),
    richDot_df = attempt(richDot(res@result, pvalue = 1, padj = 1, short = TRUE)),
    richLollipop = attempt(richLollipop(res, pvalue = 1, padj = 1, short = TRUE)),
    richCircle = attempt(richCircle(res, pvalue = 1, padj = 1, short = TRUE)),
    richGeneDot = attempt(richGeneDot(res, fc = fc, pvalue = 1, padj = 1, short = TRUE)),
    richGeneBar = attempt(richGeneBar(res, pvalue = 1, padj = 1, short = TRUE))
  )
  for (what in names(plots)) {
    expect_false(inherits(plots[[what]], "error"), info = what)
    if (inherits(plots[[what]], "error")) next
    cats <- plot_categories(plots[[what]])
    long <- unique(cats[startsWith(cats, "a")])
    expect_equal(length(long), 2, info = what)
    # short = TRUE reached the label: the eight-word labels are wrapped.
    expect_true(length(long) > 0 && all(grepl("\n", long, fixed = TRUE)), info = what)
  }
})

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
