test_that("annotation membership normalization validates and deduplicates rows", {
  annot <- data.frame(
    GeneID = c("A", "A", " ", NA, "B"),
    Term = c("T", "T", "T", "X", ""),
    Annot = paste0("name", 1:5)
  )

  normalized <- expect_warning(
    richR:::.normalize_annotation(annot, func_name = "test"),
    "3 invalid annotation row"
  )

  expect_equal(normalized$data[, 1:2], annot[1, 1:2])
  expect_equal(normalized$data$Annot, "name1")
  expect_equal(normalized$counts[c("original_rows", "valid_rows", "invalid_rows",
                                   "duplicate_rows", "unique_rows")],
               list(original_rows = 5L, valid_rows = 2L, invalid_rows = 3L,
                    duplicate_rows = 1L, unique_rows = 1L))
  expect_error(richR:::.normalize_annotation(data.frame(GeneID = "A")), "at least 2")
  list_col <- data.frame(GeneID = I(list("A")), Term = "T")
  expect_error(richR:::.normalize_annotation(list_col),
               "atomic")
  expect_error(richR:::.normalize_annotation(data.frame(GeneID = NA, Term = " ")),
               "no valid annotation memberships")

  identity_named <- data.frame(data = c("A", "B", "C"), counts = c("T", "T", "U"))
  direct <- suppressMessages(richR:::.run_ora(
    "A", identity_named, pvalue = 1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  ))
  expect_equal(direct@result$Annotated, 2)
})

test_that("normalized constructor annotations do not warn again during analysis", {
  annotation <- expect_warning(
    richR::buildAnnotFromList(list(T = c("A", NA, " "))),
    "excluded 2 invalid annotation row"
  )
  expect_no_warning(suppressMessages(richR::enrich(
    "A", annotation, pvalue = 1, minSize = 1, minGSSize = 1,
    keepRich = FALSE
  )))
})

test_that("ORA membership counts and keepRich filtering use both bounds", {
  annot <- data.frame(
    GeneID = c("A", "A", "B", "C", "D", "E", "F", "G"),
    Term = c("T1", "T1", "T1", "T2", "T2", "T3", "T4", "T5")
  )
  result <- suppressMessages(richR:::enrich_internal(
    "A", annot, pvalue = 1, minSize = 1, maxSize = 10,
    minGSSize = 1, maxGSSize = 10, keepRich = FALSE
  ))
  expect_equal(result@result$Annotated[result@result$Annot == "T1"], 2)

  fixture <- data.frame(
    Annot = c("min_only", "gs_only", "both", "rich", "size_over", "gs_over"),
    Pvalue = rep(0.01, 6), Padj = rep(0.01, 6),
    Significant = c(2, 1, 2, 1, 11, 1),
    Annotated = c(1, 2, 2, 1, 1, 11),
    RichFactor = c(0.5, 0.5, 0.5, 1, 1, 1)
  )
  false <- richR:::.filter_ora_result(fixture, 1, NULL, 2, 10, 2, 10, FALSE)
  true <- richR:::.filter_ora_result(fixture, 1, NULL, 2, 10, 2, 10, TRUE)
  expect_equal(false$Annot, "both")
  expect_equal(true$Annot, c("both", "rich"))
  expect_length(richR:::.filter_ora_result(
    transform(fixture, Pvalue = 1), 1, NULL, 2, 10, 2, 10, TRUE
  )$Annot, 0)
})

test_that("GO membership normalization precedes ontology selection", {
  annot <- data.frame(
    GeneID = c("A", "B", "C", "D", "E", "F", NA),
    GOALL = c("BP1", "BP1", "BP2", "BP2", "MF1", "MF1", "BP1"),
    Annot = c("BP term 1", "BP term 1", "BP term 2", "BP term 2",
              "MF term", "MF term", "BP term 1"),
    ONTOLOGYALL = c("BP", "BP", "BP", "BP", "MF", "MF", "BP")
  )

  result <- expect_warning(suppressMessages(richR::richGO(
    c("A", "B"), godata = annot, ontology = "BP", pvalue = 1,
    minSize = 1, minGSSize = 1, keepRich = FALSE
  )), "excluded 1 invalid annotation row")

  expect_equal(result@result$Annot, "BP1")
  expect_equal(result@result$Annotated, 2)
  expect_equal(result@result$Pvalue, stats::phyper(1, 2, 2, 2, lower.tail = FALSE))
})

test_that("ORA correction covers positive overlaps before reporting size filters", {
  annot <- data.frame(
    GeneID = c("A", "B", "A", "C", "B", "C", "A", "D", "E", "F"),
    Term = c("T1", "T1", "T2", "T2", "T3", "T3", "T4", "T4", "T5", "T5")
  )
  full <- suppressMessages(richR:::enrich_internal(
    c("A", "B", "C"), annot, pvalue = 1, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  ))
  filtered <- suppressMessages(richR:::enrich_internal(
    c("A", "B", "C"), annot, pvalue = 1, minSize = 2,
    minGSSize = 1, keepRich = FALSE
  ))

  expect_setequal(full@result$Annot, c("T1", "T2", "T3", "T4"))
  expect_equal(full@result$Padj, stats::p.adjust(full@result$Pvalue, method = "BH"))
  expect_setequal(filtered@result$Annot, c("T1", "T2", "T3"))
  expect_equal(filtered@result$Padj,
               full@result$Padj[match(filtered@result$Annot, full@result$Annot)])
})

test_that("named gene-set and batch lists reject invalid identities before processing", {
  duplicate_sets <- structure(list(c("A"), c("B")), names = c("T", "T"))
  expect_error(richR::buildAnnotFromList(duplicate_sets), "duplicate.*T")
  expect_error(richR::buildAnnotFromList(list(" " = c("A"))), "position 1")
  expect_error(richR::buildAnnotFromList(list()), "non-empty")
  expect_error(richR::buildAnnotFromList(setNames(list("A"), NA_character_)), "position 1")
  expect_error(richR::buildAnnotFromList(setNames(list("A"), "")), "position 1")
  expect_error(richR::buildAnnotFromList(list(Empty = character())), "No valid gene sets")
  mixed <- suppressMessages(richR::buildAnnotFromList(
    list(Empty = character(), Valid = c("A", "B"))
  ))
  expect_equal(unique(mixed@annot$Term), "Valid")
  expect_invisible(richR:::.validate_named_list(
    setNames(list("A", "B"), c("T", " T")), "gene_sets"
  ))

  duplicate_lists <- structure(list(c("A"), c("B")), names = c("G", "G"))
  expect_error(richR::batchEnrich(duplicate_lists, data.frame(), verbose = FALSE),
               "duplicate.*G")
  expect_error(richR::batchEnrich(list(), data.frame(), verbose = FALSE), "non-empty")
  expect_error(richR::batchEnrich(list(" " = c("A")), data.frame(), verbose = FALSE),
               "position 1")
  expect_error(richR::batchEnrich(list(c("A")), data.frame(), verbose = FALSE), "named")
  expect_error(richR::batchEnrich(setNames(list("A"), NA_character_), data.frame(),
                                  verbose = FALSE), "position 1")
})

test_that("GMT duplicate terms merge memberships once and blank terms fail with line number", {
  gmt <- tempfile(fileext = ".gmt")
  on.exit(unlink(gmt), add = TRUE)
  writeLines("single\tdesc\tA\tB", gmt)
  single <- expect_no_warning(richR::readGMT(gmt))
  expect_equal(single@annot$GeneID, c("A", "B"))

  writeLines(c("T\tdesc\tA\tB", "T\tdesc2\tB\tC", "S\tdesc\tD"), gmt)

  annot <- expect_warning(richR::readGMT(gmt), "duplicate.*T")
  expect_equal(annot@annot$GeneID[annot@annot$Term == "T"], c("A", "B", "C"))
  expect_true(all(c("S", "T") %in% annot@annot$Term))

  writeLines(c("T\tdesc\tA", " \tdesc\tB"), gmt)
  expect_error(richR::readGMT(gmt), "line 2")
})

test_that("GSEA adjusts the complete returned family with the requested method", {
  old_param <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_param), add = TRUE)

  stats <- setNames(seq(-3, 3, length.out = 100), paste0("g", 1:100))
  pathways <- list(
    low = paste0("g", 1:15),
    mid = paste0("g", 41:55),
    high = paste0("g", 86:100)
  )
  annot <- data.frame(
    GeneID = unlist(pathways, use.names = FALSE),
    Term = rep(names(pathways), lengths(pathways))
  )
  expect_error(richR:::richGSEA_internal(stats, annot, padj.method = "invalid"),
               "padj.method")

  set.seed(42)
  raw <- suppressWarnings(fgsea::fgseaMultilevel(
    pathways, stats, minSize = 5, maxSize = 500
  ))
  for (method in c("BH", "bonferroni", "none", "BY")) {
    expected <- stats::p.adjust(raw$pval, method = method)
    cutoff <- 0.01
    set.seed(42)
    result <- suppressWarnings(suppressMessages(richR:::richGSEA_internal(
      stats, annot, pvalue = 1, padj = cutoff, minSize = 5,
      padj.method = method, table = FALSE
    )))
    selected <- expected < cutoff
    expect_setequal(result@result$pathway, raw$pathway[selected])
    expect_equal(result@result$padj,
                 expected[match(result@result$pathway, raw$pathway)], tolerance = 1e-12)
    expect_equal(result@pAdjustMethod, method)
    expect_true(any(selected) && any(!selected))
  }

  set.seed(42)
  default_result <- suppressWarnings(suppressMessages(richR:::richGSEA_internal(
    stats, annot, pvalue = 1, padj = 1, minSize = 5, table = FALSE
  )))
  expect_equal(default_result@result$padj,
               stats::p.adjust(raw$pval, "BH")[match(default_result@result$pathway,
                                                      raw$pathway)], tolerance = 1e-12)
  expect_equal(default_result@pAdjustMethod, "BH")
})

test_that("GSEA adjusts backend NA p-values with the full returned family", {
  stats <- c(A = 2, B = 1, C = -1, D = -2)
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("T1", "T1", "T2", "T2"))
  raw_p <- c(0.01, NA_real_, 0.04)
  raw_pathways <- c("T1", "T2", "T3")
  result <- testthat::with_mocked_bindings(
    fgseaMultilevel = function(...) data.frame(
      pathway = raw_pathways, pval = raw_p, NES = c(1, 0, -1), size = c(2L, 2L, 1L),
      leadingEdge = I(list("A", "B", "C"))
    ),
    richR::richGSEA(stats, annot, padj = 1, minSize = 1, table = FALSE),
    .package = "richR"
  )

  adjusted <- stats::p.adjust(raw_p, method = "BH")
  expect_equal(adjusted, c(0.02, NA_real_, 0.04))
  expect_equal(result@result$pathway, c("T1", "T3"))
  expect_equal(result@result$padj, adjusted[c(1, 3)])
  expect_equal(richR::enrichmentInfo(result)$testing$tested_count, 3L)
})

test_that("batchEnrich reports empty groups and preserves successful groups", {
  annot <- data.frame(GeneID = LETTERS[1:20],
                      Term = c(rep("T", 2), rep("U", 18)))
  batch <- expect_warning(suppressMessages(richR::batchEnrich(
    list(empty = character(), successful = c("A", "B")), annot, verbose = FALSE,
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  )), "failed for 'empty'")
  expect_named(batch$results, "successful")
  expect_equal(unique(batch$combined$group), "successful")
})

test_that("batchEnrich errors when all validly named groups are empty", {
  annot <- data.frame(GeneID = LETTERS[1:20],
                      Term = c(rep("T", 2), rep("U", 18)))
  expect_warning(
    expect_error(richR::batchEnrich(
      list(empty1 = character(), empty2 = character()), annot,
      verbose = FALSE, pvalue = 1
    ), "all enrichment analyses failed"),
    "failed for 'empty1'"
  )
})

test_that("zero-variance ORA z-scores are NA", {
  stats <- richR:::.compute_ora_stats(c(T = 1), c(T = 1), N = 1, n = 1)
  expect_equal(stats$rhs, c(T = 1))
  expect_identical(stats$zscore, c(T = NA_real_))
})
