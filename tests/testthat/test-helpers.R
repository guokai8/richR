test_that(".clean.char removes newlines", {
  expect_equal(richR:::.clean.char("foo\nbar"), "foo bar")
  expect_equal(richR:::.clean.char("no newlines"), "no newlines")
})

test_that(".filter_ora_result filters by pvalue", {
  df <- data.frame(
    Annot = c("A", "B", "C"),
    Pvalue = c(0.01, 0.05, 0.1),
    Padj = c(0.03, 0.08, 0.15),
    Significant = c(5, 3, 2),
    Annotated = c(100, 50, 30),
    RichFactor = c(0.05, 0.06, 0.067)
  )
  result <- richR:::.filter_ora_result(df, pvalue = 0.06, padj = NULL,
                                        minSize = 1, maxSize = 500,
                                        minGSSize = 10, maxGSSize = 500,
                                        keepRich = TRUE)
  expect_true(nrow(result) <= 2)
  expect_true(all(result$Pvalue < 0.06))
})

test_that(".filter_ora_result filters by padj when provided", {
  df <- data.frame(
    Annot = c("A", "B", "C"),
    Pvalue = c(0.01, 0.05, 0.1),
    Padj = c(0.03, 0.08, 0.15),
    Significant = c(5, 3, 2),
    Annotated = c(100, 50, 30),
    RichFactor = c(0.05, 0.06, 0.067)
  )
  result <- richR:::.filter_ora_result(df, pvalue = 0.05, padj = 0.1,
                                        minSize = 1, maxSize = 500,
                                        minGSSize = 10, maxGSSize = 500,
                                        keepRich = TRUE)
  expect_true(all(result$Padj < 0.1))
})

test_that(".build_detail creates gene-term detail from vector input", {
  res_df <- data.frame(
    Annot = c("T1", "T2"),
    Term = c("Term1", "Term2"),
    Pvalue = c(0.01, 0.05),
    Padj = c(0.02, 0.08),
    GeneID = c("A,B", "C")
  )
  detail <- richR:::.build_detail(res_df, x = c("A", "B", "C"), sep = ",")
  expect_s3_class(detail, "data.frame")
  expect_true("GeneID" %in% colnames(detail))
  expect_true(nrow(detail) >= 3)
})

test_that("vec_to_df creates correct data.frame", {
  v <- c(a = 1, b = 2, c = 3)
  df <- vec_to_df(v, c("name", "value"))
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_equal(colnames(df), c("name", "value"))
})
