encoding_identity_pair <- function(text = "té") {
  utf8 <- enc2utf8(text)
  latin1 <- iconv(utf8, from = "UTF-8", to = "latin1")
  Encoding(latin1) <- "latin1"
  c(utf8, latin1)
}

test_that("native groups and counts use text identity across encodings", {
  pair <- encoding_identity_pair()
  expect_identical(length(unique(pair)), 1L)
  counted <- list(pair = pair, long = c(pair, "B", pair),
                  empty = character(), missing = c(NA_character_, "NA", NA),
                  other = 1:2)
  expect_identical(richR:::name_table(counted),
                   c(pair = 1L, long = 2L, empty = 0L, missing = 2L,
                     other = NA_integer_))
  grouped <- richR:::sf(data.frame(value = c("g1", "g2", "g3"),
                                  key = c(pair, "B")))
  expect_identical(grouped, setNames(list(c("g1", "g2"), "g3"),
                                    c(pair[1], "B")))
  reversed <- richR:::reverseList(list(g1 = pair[1], g2 = pair[2], g3 = "B"))
  expect_identical(reversed, setNames(list(c("g1", "g2"), "g3"),
                                     c(pair[1], "B")))
})

test_that("encoding-equivalent ORA terms match a normalized analysis", {
  pair <- encoding_identity_pair()
  annot <- data.frame(GeneID = c("g1", "g2", "g3", "g4"),
                      Term = c(pair, "B", "B"),
                      Annot = c("label", "label", "B", "B"))
  normalized <- annot
  normalized$Term <- enc2utf8(normalized$Term)
  analyze <- function(a) suppressMessages(richR::enrich(
    c("g1", "g2"), a, pvalue = 1, padj = NULL, minSize = 1, minGSSize = 1
  ))
  mixed <- analyze(annot)@result
  control <- analyze(normalized)@result
  expect_equal(mixed, control)
  expect_equal(mixed$Annotated, 2)
  expect_equal(mixed$Significant, 2)
  expect_identical(mixed$GeneID, "g1,g2")
  expect_equal(mixed$Pvalue, stats::phyper(1, 2, 2, 2, lower.tail = FALSE))
})

test_that("native keys preserve ASCII order and missing-value distinctions", {
  keys <- c("B", NA, "NA", "A", "B", NA)
  values <- paste0("g", seq_along(keys))
  expected <- setNames(list(values[c(1, 5)], values[c(2, 6)], "g3", "g4"),
                       c("B", NA, "NA", "A"))
  expect_identical(richR:::sf(data.frame(values, keys)), expected)
  expect_identical(richR:::reverseList(setNames(as.list(keys), values)), expected)
  expect_identical(richR:::reverseList(list(A = c("x", "x"), B = "x")),
                   list(x = c("A", "A", "B")))
  expect_identical(richR:::sf(data.frame(values = character(), keys = character())),
                   setNames(list(), character()))
  expect_identical(richR:::reverseList(list()), setNames(list(), character()))
  expect_identical(richR:::name_table(list(c(NA_character_, NA_character_))), 1L)
})

test_that("normalized native keys remain live during bounded GC stress", {
  # Keep assertions outside torture mode so only the native calls are stressed.
  stress <- function(pair) {
    input <- data.frame(value = c("g1", "g2"), key = pair)
    lhs <- list(g1 = pair[1], g2 = pair[2])
    old <- gctorture2(1L)
    on.exit(gctorture2(old), add = TRUE)
    list(grouped = richR:::sf(input), reversed = richR:::reverseList(lhs),
         counted = richR:::name_table(list(c(pair, pair))))
  }
  for (i in seq_len(3L)) {
    pair <- encoding_identity_pair(paste0("té", i))
    expected <- setNames(list(c("g1", "g2")), pair[1])
    result <- stress(pair)
    expect_identical(result$grouped, expected)
    expect_identical(result$reversed, expected)
    expect_identical(result$counted, 1L)
  }
})

test_that("Latin-1-only keys stay protected through output allocation", {
  latin1 <- iconv(paste0("té-protected-", seq_len(24L)),
                 from = "UTF-8", to = "latin1")
  Encoding(latin1) <- "latin1"
  input <- data.frame(value = rep(c("g1", "g2"), 4L),
                      key = rep(latin1[1:4], each = 2L))
  lhs <- setNames(as.list(rep(latin1[5:12], each = 2L)),
                  rep(c("g1", "g2"), 8L))
  counted <- lapply(latin1[13:24], rep, times = 3L)
  collect <- function() {
    old <- gctorture2(1L)
    on.exit(gctorture2(old), add = TRUE)
    list(grouped = richR:::sf(input), reversed = richR:::reverseList(lhs),
         counted = richR:::name_table(counted))
  }
  gc()
  result <- collect()
  expect_identical(result$grouped,
                   setNames(rep(list(c("g1", "g2")), 4L), enc2utf8(latin1[1:4])))
  expect_identical(result$reversed,
                   setNames(rep(list(c("g1", "g2")), 8L), enc2utf8(latin1[5:12])))
  expect_identical(result$counted, rep(1L, 12L))
})
