compare_test_annotation <- function() {
  data.frame(
    GeneID = LETTERS[1:6],
    Annot = rep(c("T1", "T2", "T3"), each = 2),
    stringsAsFactors = FALSE
  )
}

compare_test_result <- function(genes, cutoff = 0.1) {
  suppressMessages(richR::enrich(
    genes, compare_test_annotation(), pvalue = cutoff, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  ))
}

test_that("compareResult preserves default and explicit group identity", {
  first <- compare_test_result(c("A", "B"))
  second <- compare_test_result(c("C", "D"))

  unnamed <- richR::compareResult(list(first, second), pvalue = 0.1)
  expect_identical(unique(unnamed$group), c("Group 1", "Group 2"))
  expect_identical(unnamed$Annot, c("T1", "T2"))
  expect_error(richR::compareResult(setNames(list(first, second), c("same", "same"))), "same")
  expect_error(richR::compareResult(setNames(list(first, second), c("ok", NA))), "missing|position")
  expect_error(richR::compareResult(setNames(list(first, second), c("ok", "  "))), "blank|position")
  expect_error(richR::compareResult(list(named = first, second)), "missing or blank names.*position 2")
  expect_error(richR::compareResult(list()), "non-empty")
  expect_error(richR::compareResult("not a list"), "list")
  expect_error(richR::compareResult(list(bad = 42)), "bad.*richResult.*data.frame.*NULL")
})

test_that("compareResult keeps strict raw and adjusted p-value boundaries", {
  boundary <- data.frame(
    Annot = c("p-boundary", "padj-boundary", "inside", "na-p", "na-padj"),
    Term = LETTERS[1:5], Annotated = rep(2, 5), Significant = rep(1, 5),
    RichFactor = rep(0.5, 5), FoldEnrichment = rep(1, 5), zscore = rep(0, 5),
    Pvalue = c(0.05, 0.01, 0.049, NA, 0.01),
    Padj = c(0.01, 1, 0.99, 0.01, NA),
    GeneID = LETTERS[1:5], stringsAsFactors = FALSE
  )

  pfiltered <- richR::compareResult(list(g = boundary), pvalue = 0.05)
  expect_identical(pfiltered$Annot, "inside")
  padj_filtered <- richR::compareResult(list(g = boundary), padj = 1)
  expect_identical(padj_filtered$Annot, c("p-boundary", "inside"))
})

test_that("compareResult reports missing required columns with the group name", {
  incomplete <- data.frame(Annot = "T1", Pvalue = 1 / 15)
  expect_error(
    richR::compareResult(list(incomplete_group = incomplete)),
    "group 'incomplete_group'.*required result column.*GeneID"
  )
})

test_that("compareResult preserves real extra columns for mixed result inputs", {
  first <- compare_test_result(c("A", "B"))
  second <- compare_test_result(c("C", "D"))
  first@result$extra_first <- "first"
  second_frame <- second@result
  second_frame$extra_second <- 7L

  compared <- richR::compareResult(
    list(first = first, second = second_frame), pvalue = 0.1
  )

  expect_identical(as.character(compared$group), c("first", "second"))
  expect_identical(compared$Annot, c("T1", "T2"))
  expect_identical(compared$extra_first, c("first", NA_character_))
  expect_identical(compared$extra_second, c(NA_integer_, 7L))
})

test_that("include.all preserves axes without relaxing the default cutoff", {
  first <- compare_test_result(c("A", "B"))
  result <- richR::compareResult(list(first = first), include.all = TRUE)

  expect_equal(first@result$Pvalue, 1 / 15, tolerance = 1e-12)
  expect_equal(nrow(result), 0L)
  expect_identical(levels(result$group), "first")
})

test_that("include.all preserves axes without manufacturing result rows", {
  first <- compare_test_result(c("A", "B"))
  empty <- compare_test_result(c("A", "C"))
  compared <- richR::compareResult(
    list(first = first, empty = empty, missing = NULL), pvalue = 0.1,
    include.all = TRUE
  )

  expect_identical(levels(compared$group), c("first", "empty", "missing"))
  expect_identical(unique(as.character(compared$group)), "first")
  expect_equal(nrow(compared), 1L)
  expect_false(any(compared$Pvalue == 1))

  plot <- richR::richCompareDot(compared, include.all = TRUE)
  built <- ggplot2::ggplot_build(plot)
  expect_identical(built$layout$panel_params[[1]]$x$get_labels(), c("first", "empty", "missing"))
  expect_equal(nrow(built$data[[1]]), 1L)

  filtered <- richR::richCompareDot(compared, include.all = TRUE, pathway = "T1", level = "Term")
  filtered_build <- ggplot2::ggplot_build(filtered)
  expect_identical(filtered_build$layout$panel_params[[1]]$x$get_labels(), c("first", "empty", "missing"))
  expect_equal(nrow(filtered_build$data[[1]]), 1L)

  dropped <- richR::richCompareDot(compared, pvalue = 0.1, pathway = "T1", level = "Term")
  dropped_build <- ggplot2::ggplot_build(dropped)
  expect_identical(dropped_build$layout$panel_params[[1]]$x$get_labels(), "first")
  expect_equal(nrow(dropped_build$data[[1]]), 1L)
})

test_that("all-empty and all-NULL comparisons keep typed frames and axes", {
  empty <- compare_test_result(c("A", "C"))
  all_empty <- richR::compareResult(list(empty = empty, missing = NULL), include.all = TRUE)
  expect_equal(nrow(all_empty), 0L)
  expect_identical(levels(all_empty$group), c("empty", "missing"))
  expect_type(all_empty$Pvalue, "double")
  expect_type(all_empty$Annot, "character")
  empty_plot <- richR::richCompareDot(all_empty, include.all = TRUE)
  empty_build <- ggplot2::ggplot_build(empty_plot)
  expect_identical(empty_build$layout$panel_params[[1]]$x$get_labels(), c("empty", "missing"))
  expect_equal(nrow(empty_build$data[[1]]), 0L)
  filtered_empty_build <- ggplot2::ggplot_build(richR::richCompareDot(
    all_empty, include.all = TRUE, pathway = "absent", level = "Term"
  ))
  expect_identical(filtered_empty_build$layout$panel_params[[1]]$x$get_labels(), c("empty", "missing"))
  expect_equal(nrow(filtered_empty_build$data[[1]]), 0L)
  empty_pdf <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(empty_pdf, plot = empty_plot, width = 5, height = 4)
  expect_true(file.exists(empty_pdf))
  expect_gt(file.info(empty_pdf)$size, 0)
  unlink(empty_pdf)

  null_comparison <- richR::compareResult(list(first = NULL, second = NULL), include.all = TRUE)
  expect_equal(nrow(null_comparison), 0L)
  expect_identical(levels(null_comparison$group), c("first", "second"))
  expect_type(null_comparison$Pvalue, "double")
  null_plot <- richR::richCompareDot(null_comparison, include.all = TRUE)
  expect_identical(ggplot2::ggplot_build(null_plot)$layout$panel_params[[1]]$x$get_labels(), c("first", "second"))

  wide_names <- paste0("group_", seq_len(12))
  wide <- richR::compareResult(
    setNames(rep(list(empty), length(wide_names)), wide_names), include.all = TRUE
  )
  wide_plot <- ggplot2::ggplot_build(richR::richCompareDot(wide, include.all = TRUE))
  expect_identical(wide_plot$layout$panel_params[[1]]$x$get_labels(), wide_names)
  expect_equal(nrow(wide_plot$data[[1]]), 0L)
})
