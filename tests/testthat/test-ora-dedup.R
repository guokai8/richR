test_that("ORA treats repeated vector genes as one input gene", {
  annot <- data.frame(
    GeneID = c("A", "B", "C", "D", "E", "F"),
    Term = c("T1", "T1", "T2", "T2", "T3", "T3")
  )

  unique_result <- suppressMessages(richR::enrich(
    c("A", "B"), annot, pvalue = 1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))
  duplicate_result <- suppressWarnings(suppressMessages(richR::enrich(
    c("A", "B", "A"), annot, pvalue = 1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  )))

  expect_equal(duplicate_result@result, unique_result@result)
  expect_equal(duplicate_result@gene, c("A", "B"))
  expect_equal(duplicate_result@genenumber, 2)
  expect_equal(duplicate_result@detail, unique_result@detail)

  annot_object <- new("Annot", species = "custom", anntype = "custom",
                      keytype = "SYMBOL", annot = annot)
  object_result <- suppressWarnings(suppressMessages(richR::enrich(
    c("A", "B", "A"), annot_object, pvalue = 1, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  )))
  expect_equal(object_result@result, unique_result@result)
  expect_equal(object_result@gene, c("A", "B"))
  expect_equal(object_result@genenumber, 2)
})

test_that("duplicated ORA input preserves first-seen order and the exact tail probability", {
  annot <- data.frame(
    GeneID = c("A", "B", "C", "D", "E", "F"),
    Term = c("T1", "T1", "T2", "T2", "T3", "T3")
  )
  result <- suppressWarnings(suppressMessages(richR::enrich(
    c("B", "A", "B"), annot, pvalue = 1, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  )))

  expect_identical(result@gene, c("B", "A"))
  expect_identical(result@genenumber, 2L)
  expect_equal(result@result$Pvalue[result@result$Annot == "T1"],
               1 / 15, tolerance = 1e-15)
  expect_equal(result@result$Pvalue[result@result$Annot == "T1"],
               stats::phyper(1, 2, 4, 2, lower.tail = FALSE), tolerance = 1e-15)
})

test_that("ORA preserves data-frame input detail metadata", {
  annot <- data.frame(
    GeneID = c("A", "B", "C", "D", "E", "F"),
    Term = c("T1", "T1", "T2", "T2", "T3", "T3")
  )
  input <- data.frame(logFC = c(1.5, -2), row.names = c("A", "B"))

  result <- suppressMessages(richR::enrich(
    input, annot, pvalue = 1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))

  expect_equal(result@gene, c("A", "B"))
  expect_equal(result@genenumber, 2)
  expect_equal(result@detail$logFC[match(c("A", "B"), result@detail$GeneID)],
               c(1.5, -2))
})
