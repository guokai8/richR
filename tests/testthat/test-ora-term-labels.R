test_that("KEGG falls back to term IDs for missing or blank labels", {
  annotation <- data.frame(
    GeneID = c("A", "B", "C", "D"),
    PATH = c("99999", "99999", "10000", "10000")
  )
  object <- new("Annot", species = "custom", anntype = "KEGG",
                keytype = "SYMBOL", annot = annotation)
  lookup_with_blank <- data.frame(annotation = "   ", row.names = "10000")

  run <- function(input, labels) {
    testthat::with_mocked_bindings(
      .get_kg_dat = function(builtin = TRUE) labels,
      richR::richKEGG(c("A", "C"), input, pvalue = 1, minSize = 1,
                      minGSSize = 1, keepRich = FALSE),
      .package = "richR"
    )
  }

  missing_and_blank <- run(annotation, lookup_with_blank)
  missing_and_blank_annot <- run(object, lookup_with_blank)
  for (result in list(missing_and_blank, missing_and_blank_annot)) {
    expect_setequal(result@result$Term, c("99999", "10000"))
    expect_equal(result@result$Annotated, c(2, 2))
    expect_equal(result@result$Significant, c(1, 1))
    expect_equal(richR::enrichmentInfo(result)$annotation$effective_rows, 4L)
    expect_equal(richR::enrichmentInfo(result)$testing$tested_count, 2L)
  }
  expect_equal(missing_and_blank@result$Pvalue,
               missing_and_blank_annot@result$Pvalue)

  lookup_with_label <- data.frame(annotation = "Known pathway", row.names = "10000")
  valid_lookup <- run(annotation, lookup_with_label)
  valid_lookup_annot <- run(object, lookup_with_label)
  expect_equal(valid_lookup@result$Term[valid_lookup@result$Annot == "10000"],
               "Known pathway")
  expect_equal(valid_lookup_annot@result$Term[valid_lookup_annot@result$Annot == "10000"],
               "Known pathway")
  expect_equal(valid_lookup@result$Pvalue, missing_and_blank@result$Pvalue)
  expect_equal(valid_lookup@result$Annotated, missing_and_blank@result$Annotated)
})
