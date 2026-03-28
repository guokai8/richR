context("S4 class construction and show methods")

test_that("richResult can be constructed and displayed", {
  res <- new("richResult",
    result = data.frame(
      Annot = "GO:0001", Term = "test term", Annotated = 100,
      Significant = 10, RichFactor = 0.1, FoldEnrichment = 2.0,
      zscore = 3.0, Pvalue = 0.001, Padj = 0.01,
      GeneID = "A,B,C", stringsAsFactors = FALSE
    ),
    detail = data.frame(TERM = "GO:0001", Annot = "test term",
                        GeneID = "A", Pvalue = 0.001, Padj = 0.01),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    padjCutoff = numeric(),
    genenumber = 3L,
    organism = "human",
    ontology = "BP",
    gene = c("A", "B", "C"),
    keytype = "SYMBOL",
    sep = ","
  )

  expect_s4_class(res, "richResult")
  expect_equal(nrow(res@result), 1)
  expect_output(show(res), "richResult object")
  expect_output(show(res), "Significant:.*1 terms")
})

test_that("GSEAResult can be constructed and displayed", {
  res <- new("GSEAResult",
    result = data.frame(
      pathway = "MAPK signaling", pval = 0.001, padj = 0.01,
      log2err = 0.5, ES = 0.6, NES = 1.8, size = 50,
      leadingEdge = "A,B,C", stringsAsFactors = FALSE
    ),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    padjCutoff = numeric(),
    genenumber = 100L,
    organism = "human",
    gene = c("A", "B", "C"),
    input = c(A = 2.1, B = -1.5, C = 0.3),
    ontology = "KEGG",
    keytype = "SYMBOL",
    sep = ","
  )

  expect_s4_class(res, "GSEAResult")
  expect_output(show(res), "GSEAResult object")
  expect_output(show(res), "MAPK signaling")
})

test_that("Annot can be constructed and displayed", {
  annot <- new("Annot",
    species = "human",
    anntype = "KEGG",
    keytype = "SYMBOL",
    annot = data.frame(GeneID = c("A", "B", "C"),
                       Term = c("path1", "path1", "path2"))
  )

  expect_s4_class(annot, "Annot")
  expect_output(show(annot), "Annot object")
  expect_output(show(annot), "Genes:.*3")
  expect_output(show(annot), "Terms:.*2")
})

test_that("richResult with zero results displays correctly", {
  res <- new("richResult",
    result = data.frame(Annot = character(), Term = character(),
                        Annotated = integer(), Significant = integer(),
                        RichFactor = numeric(), FoldEnrichment = numeric(),
                        zscore = numeric(), Pvalue = numeric(), Padj = numeric(),
                        GeneID = character(), stringsAsFactors = FALSE),
    detail = data.frame(),
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    padjCutoff = numeric(),
    genenumber = 100L,
    organism = "human",
    ontology = "BP",
    gene = c("X", "Y"),
    keytype = "SYMBOL",
    sep = ","
  )

  expect_output(show(res), "Significant:.*0 terms")
})
