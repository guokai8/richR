test_that("ORA uses the supplied effective universe for its hypergeometric test", {
  annot <- data.frame(
    GeneID = c("A", "B", "B", "C", "D", "E", "F"),
    Term = c("T1", "T1", "T2", "T2", "T3", "T3", "T3")
  )
  query <- c("A", "B")
  expected_p <- c(T1 = stats::phyper(1, 2, 2, 2, lower.tail = FALSE),
                  T2 = stats::phyper(0, 2, 2, 2, lower.tail = FALSE))

  res <- suppressMessages(richR::enrich(
    query, annot, universe = c("A", "B", "C", "D"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))

  expect_equal(res@result$Pvalue[match(names(expected_p), res@result$Annot)],
               unname(expected_p))
  expect_equal(res@result$Annotated[match(names(expected_p), res@result$Annot)],
               c(2, 2))
  expect_equal(res@result$Significant[match(names(expected_p), res@result$Annot)],
               c(2, 1))
  expect_equal(res@gene, query)
  expect_equal(res@genenumber, 2)
  object <- new("Annot", species = "custom", anntype = "custom",
                keytype = "SYMBOL", annot = annot)
  object_result <- suppressMessages(richR::enrich(
    query, object, universe = c("A", "B", "C", "D"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(object_result@result, res@result)
  object_sep <- suppressMessages(richR::enrich(
    query, object, universe = c("A", "B", "C", "D"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE,
    sep = "|", builtin = FALSE
  ))
  expect_equal(object_sep@result$GeneID[match("T1", object_sep@result$Annot)], "A|B")

  default <- suppressMessages(richR::enrich(
    query, annot, pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  all_genes <- suppressMessages(richR::enrich(
    query, annot, universe = c("A", "B", "C", "D", "E", "F"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(all_genes@result, default@result)
  positional <- suppressMessages(richR::enrich(
    query, annot, "", 1, NULL, NULL, "", NULL, 1, 500, 1, 500,
    FALSE, "BH", ","
  ))
  expect_equal(positional@result, default@result)
})

test_that("ORA cleans explicit universes and errors when no effective query remains", {
  annot <- data.frame(GeneID = c("A", "B", "C"), Term = c("T", "T", "U"))
  res <- suppressMessages(richR::enrich(
    c("A", "B", "X"), annot, universe = c("A", "C", "Y", "Y", NA, ""),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(res@result$Pvalue[match("T", res@result$Annot)],
               stats::phyper(0, 1, 1, 1, lower.tail = FALSE))
  expect_equal(res@gene, c("A", "B", "X"))
  expect_equal(res@genenumber, 3)

  expect_error(richR::enrich("A", annot, universe = numeric(1)),
               "atomic character")
  expect_error(richR::enrich("A", annot, universe = matrix("A")),
               "dimensions")
  expect_error(richR::enrich("A", annot, universe = character()),
               "empty")
  expect_error(richR::enrich("A", annot, universe = c(NA_character_, " ")),
               "no valid genes")
  expect_error(richR::enrich("A", annot, universe = "Z"),
               "no annotated genes")
  expect_error(richR::enrich("B", annot, universe = "A"),
               "no query genes")
})

test_that("richGO restricts the universe after selecting its ontology", {
  annot <- data.frame(
    GeneID = c("A", "B", "B", "C", "D", "E", "F"),
    GOALL = c("T1", "T1", "T2", "T2", "T3", "T3", "T3"),
    ONTOLOGYALL = c("BP", "BP", "BP", "BP", "BP", "CC", "CC"),
    Annot = c("one", "one", "two", "two", "three", "other", "other")
  )
  res <- suppressMessages(richR::richGO(
    c("A", "B"), annot, ontology = "BP", universe = c("A", "B", "C", "D", "E"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(res@result$Pvalue[match(c("T1", "T2"), res@result$Annot)],
               c(stats::phyper(1, 2, 2, 2, lower.tail = FALSE),
                 stats::phyper(0, 2, 2, 2, lower.tail = FALSE)))
  object <- new("Annot", species = "custom", anntype = "GO",
                keytype = "SYMBOL", annot = annot)
  object_result <- suppressMessages(richR::richGO(
    c("A", "B"), object, ontology = "BP", universe = c("A", "B", "C", "D", "E"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(object_result@result, res@result)
})

test_that("richKEGG applies the effective universe for data.frame and Annot inputs", {
  terms <- rownames(richR:::.get_kg_dat())[1:2]
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = rep(terms, each = 2))
  expected <- stats::phyper(1, 2, 2, 2, lower.tail = FALSE)
  res <- suppressMessages(richR::richKEGG(
    c("A", "B"), annot, universe = c("A", "B", "C", "D"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(res@result$Pvalue[match(terms[1], res@result$Annot)], expected)

  object <- new("Annot", species = "custom", anntype = "KEGG",
                keytype = "SYMBOL", annot = annot)
  object_result <- suppressMessages(richR::richKEGG(
    c("A", "B"), object, universe = c("A", "B", "C", "D"),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_equal(object_result@result, res@result)
})

test_that("batchEnrich forwards the universe to ORA methods", {
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("T1", "T1", "T2", "T2"))
  batch <- suppressMessages(richR::batchEnrich(
    list(sample = c("A", "B")), annot, fun = "enrich", verbose = FALSE,
    universe = c("A", "B", "C"), pvalue = 1, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  ))
  result <- batch$results[["sample"]]
  expect_equal(result@result$Pvalue[match("T1", result@result$Annot)],
               stats::phyper(1, 2, 1, 2, lower.tail = FALSE))
  expect_equal(result@gene, c("A", "B"))
})
