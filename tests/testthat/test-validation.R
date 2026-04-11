test_that(".validate_enrichment_input rejects empty gene vector", {
  expect_error(
    richR:::.validate_enrichment_input(character(0)),
    "empty"
  )
})

test_that(".validate_enrichment_input rejects all-NA genes", {
  expect_error(
    richR:::.validate_enrichment_input(c(NA, NA, NA)),
    "NA"
  )
})

test_that(".validate_enrichment_input warns on duplicated genes", {
  expect_warning(
    richR:::.validate_enrichment_input(c("A", "B", "A")),
    "[Dd]uplicated"
  )
})

test_that(".validate_enrichment_input rejects bad pvalue", {
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), pvalue = -1),
    "pvalue"
  )
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), pvalue = 2),
    "pvalue"
  )
})

test_that(".validate_enrichment_input rejects bad padj", {
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), padj = -1),
    "padj"
  )
})

test_that(".validate_enrichment_input rejects bad padj.method", {
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), padj.method = "foobar"),
    "padj.method"
  )
})

test_that(".validate_enrichment_input rejects minSize > maxSize", {
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), minSize = 100, maxSize = 10),
    "minSize"
  )
})

test_that(".validate_enrichment_input rejects bad annotation", {
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), annot = data.frame(x = 1)),
    "2 columns"
  )
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"), annot = data.frame()),
    "2 columns"
  )
  expect_error(
    richR:::.validate_enrichment_input(c("A", "B"),
      annot = data.frame(GeneID = character(0), Term = character(0))),
    "empty"
  )
})

test_that(".validate_enrichment_input accepts valid inputs", {
  expect_true(
    richR:::.validate_enrichment_input(c("GeneA", "GeneB", "GeneC"))
  )
  annot <- data.frame(GeneID = c("A", "B"), Term = c("T1", "T2"))
  expect_true(
    richR:::.validate_enrichment_input(c("A", "B"), annot = annot)
  )
})
