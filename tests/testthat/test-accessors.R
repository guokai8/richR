context("Accessor methods")

# Helper to create a test richResult
make_test_result <- function() {
  new("richResult",
    result = data.frame(
      Annot = c("GO:0001", "GO:0002"),
      Term = c("apoptosis", "cell cycle"),
      Annotated = c(100, 200),
      Significant = c(10, 20),
      RichFactor = c(0.1, 0.1),
      FoldEnrichment = c(2.0, 2.0),
      zscore = c(3.0, 3.0),
      Pvalue = c(0.001, 0.01),
      Padj = c(0.01, 0.05),
      GeneID = c("A,B,C", "D,E,F,G"),
      stringsAsFactors = FALSE
    ),
    detail = data.frame(),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    padjCutoff = numeric(),
    genenumber = 10L,
    organism = "human",
    ontology = "BP",
    gene = LETTERS[1:10],
    keytype = "SYMBOL",
    sep = ","
  )
}

test_that("result() returns data.frame", {
  res <- make_test_result()
  df <- result(res)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
})

test_that("head/tail work on richResult", {
  res <- make_test_result()
  expect_output(h <- head(res, 1), "Total significant")
  expect_equal(nrow(h), 1)
})

test_that("dim works on richResult", {
  res <- make_test_result()
  expect_equal(dim(res), c(2, 10))
})

test_that("$ accessor works on richResult", {
  res <- make_test_result()
  expect_equal(res$Pvalue, c(0.001, 0.01))
  expect_equal(res$Term, c("apoptosis", "cell cycle"))
})

test_that("[ subsetting works on richResult", {
  res <- make_test_result()
  expect_equal(res[1, "Term"], "apoptosis")
})

test_that("names() works on richResult", {
  res <- make_test_result()
  expect_true("Pvalue" %in% names(res))
})

test_that("as.data.frame works on richResult", {
  res <- make_test_result()
  df <- as.data.frame(res)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
})

test_that("getGenes extracts genes from a term", {
  res <- make_test_result()
  genes <- getGenes(res, "GO:0001")
  expect_equal(sort(genes), c("A", "B", "C"))
})

test_that("getGenes returns all genes when no term specified", {
  res <- make_test_result()
  genes <- getGenes(res)
  expect_equal(sort(genes), sort(LETTERS[1:10]))
})
