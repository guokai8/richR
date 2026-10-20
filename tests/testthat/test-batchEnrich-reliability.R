batch_test_annotation <- function() {
  data.frame(
    GeneID = LETTERS[1:6],
    Annot = rep(c("T1", "T2", "T3"), each = 2),
    stringsAsFactors = FALSE
  )
}

test_that("batchEnrich reports success, valid empty, and controlled failure", {
  annot <- batch_test_annotation()
  base_enrich <- richR::enrich
  groups <- list(
    success = c("A", "B"),
    empty = c("A", "C"),
    failed = "CONTROLLED_FAILURE"
  )
  result <- testthat::with_mocked_bindings(
    enrich = function(x, object, ...) {
      if (identical(x, "CONTROLLED_FAILURE")) stop("controlled failure")
      base_enrich(x, object, ...)
    },
    expect_warning(suppressMessages(richR::batchEnrich(
      groups, annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
      keepRich = FALSE, verbose = FALSE
    )), "failed for 'failed'"),
    .package = "richR"
  )

  expect_named(result, c("results", "combined", "status"))
  expect_named(result$results, c("success", "empty"))
  expect_equal(nrow(result$results$success@result), 1L)
  expect_identical(result$results$success@result$Annot, "T1")
  expect_equal(nrow(result$results$empty@result), 0L)
  expect_equal(nrow(result$combined), 1L)
  expect_identical(levels(result$combined$group), c("success", "empty"))
  expect_identical(result$status$group, c("success", "empty", "failed"))
  expect_identical(result$status$status, c("success", "empty", "failed"))
  expect_type(result$status$group, "character")
  expect_type(result$status$status, "character")
  expect_type(result$status$n_terms, "integer")
  expect_type(result$status$error, "character")
  expect_identical(result$status$n_terms, c(1L, 0L, NA_integer_))
  expect_identical(result$status$error, c(NA_character_, NA_character_, "controlled failure"))
})

test_that("all-valid-empty batches return typed combined data and group levels", {
  annot <- batch_test_annotation()
  result <- suppressMessages(richR::batchEnrich(
    list(first = c("A", "C"), second = c("A", "E")), annot,
    pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE, verbose = FALSE
  ))

  expect_named(result$results, c("first", "second"))
  expect_true(all(vapply(result$results, function(x) nrow(x@result) == 0L, logical(1))))
  expect_equal(nrow(result$combined), 0L)
  expect_identical(levels(result$combined$group), c("first", "second"))
  expect_identical(result$status$status, c("empty", "empty"))
  expect_identical(result$status$n_terms, c(0L, 0L))
  expect_identical(result$status$error, c(NA_character_, NA_character_))
  expect_type(result$combined$Annot, "character")
  expect_type(result$combined$Pvalue, "double")
})

test_that("all-failed batches retain warnings and the established error", {
  annot <- batch_test_annotation()
  testthat::with_mocked_bindings(
    enrich = function(x, object, ...) stop("controlled failure"),
    expect_warning(expect_error(richR::batchEnrich(
      list(first = "A", second = "B"), annot, verbose = FALSE
    ), "all enrichment analyses failed"), "failed for 'first'"),
    .package = "richR"
  )
})

test_that("batchEnrich rejects missing, blank, and duplicate group names", {
  annot <- batch_test_annotation()
  expect_error(richR::batchEnrich(setNames(list("A", "B"), c("first", NA)), annot), "missing|position")
  expect_error(richR::batchEnrich(setNames(list("A", "B"), c("first", "  ")), annot), "blank|position")
  expect_error(richR::batchEnrich(list(same = "A", same = "B"), annot), "same")
  expect_error(richR::batchEnrich(list(named = "A", "B"), annot), "missing or blank names.*position 2")
})
