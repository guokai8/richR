context("Input validation")

test_that("empty gene input raises error", {
  expect_error(
    richR:::.validateGeneInput(character(0), func_name = "test"),
    "input gene list is empty"
  )
  expect_error(
    richR:::.validateGeneInput(c(NA, NA), func_name = "test"),
    "input gene list is empty"
  )
})

test_that("NAs and duplicates are removed with message", {
  expect_message(
    result <- richR:::.validateGeneInput(c("A", "B", NA, "A"), func_name = "test"),
    "removed"
  )
  expect_equal(result, c("A", "B"))
})

test_that("zero-overlap with annotation raises error", {
  annot <- data.frame(GeneID = c("X", "Y", "Z"), Term = c("t1", "t1", "t2"))
  expect_error(
    richR:::.validateGeneInput(c("A", "B"), annotation = annot, func_name = "test"),
    "none of the input genes were found"
  )
})

test_that("partial overlap reports mapping rate", {
  annot <- data.frame(GeneID = c("A", "X", "Y"), Term = c("t1", "t1", "t2"))
  expect_message(
    richR:::.validateGeneInput(c("A", "B", "C"), annotation = annot, func_name = "test"),
    "1/3"
  )
})

test_that("parameter validation catches invalid ranges", {
  expect_error(
    richR:::.validateParams(pvalue = -0.1, func_name = "test"),
    "pvalue must be between"
  )
  expect_error(
    richR:::.validateParams(minSize = 100, maxSize = 10, func_name = "test"),
    "minSize.*cannot be greater"
  )
  expect_error(
    richR:::.validateParams(minGSSize = 500, maxGSSize = 10, func_name = "test"),
    "minGSSize.*cannot be greater"
  )
})

test_that("empty result handler provides informative message", {
  empty_df <- data.frame(Annot = character(), Pvalue = numeric())
  expect_message(
    richR:::.handleEmptyResult(empty_df, func_name = "test"),
    "no significant terms found"
  )
})

test_that("GSEA input validation catches non-numeric input", {
  expect_error(
    richR:::richGSEA_internal(c("A", "B"), data.frame(GeneID = "A", Term = "t1")),
    "named numeric vector"
  )
})

test_that("GSEA input validation handles NA and Inf", {
  x <- c(A = 1.5, B = NA, C = Inf, D = -0.5)
  # Should not error -- NAs and Infs are silently removed
  # We can't call the full richGSEA_internal without fgsea, but we can test the validation path
  expect_error(
    richR:::richGSEA_internal(c(A = NA, B = Inf), data.frame(GeneID = "A", Term = "t1")),
    "no valid"
  )
})
