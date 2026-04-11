test_that("richResult class can be instantiated", {
  res_df <- data.frame(Annot = "GO:001", Term = "test", Annotated = 100,
                       Significant = 10, Pvalue = 0.01, Padj = 0.05,
                       GeneID = "A,B,C")
  obj <- new("richResult",
    result = res_df,
    detail = data.frame(),
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
  expect_s4_class(obj, "richResult")
  expect_equal(obj@organism, "human")
  expect_equal(obj@sep, ",")
})

test_that("GSEAResult class can be instantiated", {
  res_df <- data.frame(pathway = "path1", pval = 0.01, padj = 0.05,
                       ES = 0.5, NES = 1.2, size = 50,
                       leadingEdge = "A,B", stringsAsFactors = FALSE)
  obj <- new("GSEAResult",
    result = res_df,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    padjCutoff = numeric(),
    genenumber = 100L,
    organism = "human",
    gene = c("A", "B", "C"),
    input = c(A = 1.5, B = -0.3, C = 0.8),
    ontology = "KEGG",
    keytype = "SYMBOL",
    sep = ","
  )
  expect_s4_class(obj, "GSEAResult")
  expect_equal(obj@ontology, "KEGG")
})

test_that("Annot class can be instantiated", {
  annot_df <- data.frame(GeneID = c("A", "B"), Term = c("T1", "T2"),
                         Annot = c("GO:001", "GO:002"))
  obj <- new("Annot",
    species = "human",
    anntype = "GO",
    keytype = "SYMBOL",
    annot = annot_df
  )
  expect_s4_class(obj, "Annot")
  expect_equal(obj@species, "human")
})

test_that(".make_richResult creates valid object", {
  res_df <- data.frame(Annot = "GO:001", Term = "test", Annotated = 100,
                       Significant = 10, Pvalue = 0.01, Padj = 0.05,
                       GeneID = "A,B,C")
  detail <- data.frame(TERM = "GO:001", Annot = "test", GeneID = "A",
                       Pvalue = 0.01, Padj = 0.05)
  obj <- richR:::.make_richResult(res_df, detail, pvalue = 0.05, padj = NULL,
                                   padj.method = "BH", input = c("A", "B", "C"),
                                   organism = "human", ontology = "BP",
                                   keytype = "SYMBOL", sep = ",")
  expect_s4_class(obj, "richResult")
  expect_equal(obj@genenumber, 3L)
  expect_equal(obj@padjCutoff, numeric())
})

test_that("as.data.frame works for richResult", {
  res_df <- data.frame(Annot = "GO:001", Term = "test", Annotated = 100,
                       Significant = 10, Pvalue = 0.01, Padj = 0.05,
                       GeneID = "A,B,C")
  obj <- new("richResult",
    result = res_df, detail = data.frame(),
    pvalueCutoff = 0.05, pAdjustMethod = "BH",
    padjCutoff = numeric(), genenumber = 3L,
    organism = "human", ontology = "BP",
    gene = c("A", "B", "C"), keytype = "SYMBOL", sep = ","
  )
  df <- as.data.frame(obj)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_equal(df$Annot, "GO:001")
})
