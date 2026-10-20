batch_binding_annotation <- function() {
  data.frame(
    GeneID = LETTERS[1:6],
    Annot = rep(c("T1", "T2", "T3"), each = 2),
    stringsAsFactors = FALSE
  )
}

test_that("batchEnrich binds exact rows and extra columns from real ORA results", {
  annot <- batch_binding_annotation()
  first <- suppressMessages(richR::enrich(
    c("A", "B"), annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))
  second <- suppressMessages(richR::enrich(
    c("C", "D"), annot, pvalue = 0.1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))

  expect_identical(first@result$Annot, "T1")
  expect_identical(second@result$Annot, "T2")
  expect_equal(first@result$Pvalue, 1 / 15, tolerance = 1e-12)
  expect_equal(second@result$Pvalue, 1 / 15, tolerance = 1e-12)
  expect_equal(nrow(first@result), 1L)
  expect_equal(nrow(second@result), 1L)

  first@result$extra_first <- "first"
  second@result$extra_second <- 7L
  batch <- testthat::with_mocked_bindings(
    enrich = function(x, object, ...) {
      if (identical(x, c("A", "B"))) first else second
    },
    suppressMessages(richR::batchEnrich(
      list(first = c("A", "B"), second = c("C", "D")), annot,
      pvalue = 0.1, minSize = 1, minGSSize = 1,
      keepRich = FALSE, verbose = FALSE
    )),
    .package = "richR"
  )

  expect_named(batch$results, c("first", "second"))
  expect_equal(nrow(batch$combined), 2L)
  expect_equal(batch$combined$Pvalue, c(first@result$Pvalue, second@result$Pvalue))
  expect_identical(batch$combined$Annot, c("T1", "T2"))
  expect_identical(batch$combined$extra_first, c("first", NA_character_))
  expect_identical(batch$combined$extra_second, c(NA_integer_, 7L))
  expect_equal(nrow(richR::compareResult(batch$results)), 0L)
})
