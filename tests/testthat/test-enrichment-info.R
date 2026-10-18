legacy_rich_result <- function() {
  new("richResult", result = data.frame(), detail = data.frame(),
      pvalueCutoff = 0.05, pAdjustMethod = "BH", padjCutoff = numeric(),
      genenumber = 1, organism = character(), ontology = character(),
      gene = "A", keytype = "SYMBOL", sep = ",")
}

legacy_gsea_result <- function() {
  new("GSEAResult", result = data.frame(), pvalueCutoff = 0.05,
      pAdjustMethod = "BH", padjCutoff = numeric(), genenumber = 1,
      organism = character(), gene = "A", input = c(A = 1),
      ontology = character(), keytype = "SYMBOL", sep = ",")
}

test_that("enrichmentInfo reports explicit-background input and membership counts", {
  annot <- data.frame(
    GeneID = c("A", "A", "B", NA, ""),
    Term = c("T", "T", "T", "T", "T")
  )
  result <- suppressWarnings(suppressMessages(richR::enrich(
    c("A", "A", "B", "X", NA, ""), annot, universe = "A",
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  )))
  info <- richR::enrichmentInfo(result)

  expect_equal(info$schema_version, 1L)
  expect_true(info$available)
  expect_equal(info$analysis, "ORA")
  expect_equal(unname(unlist(info$input[c("original_count", "valid_count", "unique_count",
                                          "duplicate_count", "invalid_count")])),
               c(6L, 4L, 3L, 1L, 2L))
  expect_equal(info$input$ids, c("A", "B", "X"))
  expect_equal(info$input$mapped_ids, c("A", "B"))
  expect_equal(info$input$unmapped_ids, "X")
  expect_equal(info$input$effective_ids, "A")
  expect_equal(info$input$excluded_by_universe_ids, "B")
  expect_equal(unname(unlist(info$annotation[c("original_rows", "valid_rows", "invalid_rows",
                                               "duplicate_rows", "unique_rows", "eligible_rows",
                                               "effective_rows")])),
               c(5L, 3L, 2L, 1L, 2L, 2L, 1L))
  expect_equal(info$background$mode, "explicit")
  expect_equal(info$background$effective_count, 1L)
  expect_equal(info$background$effective_ids, "A")
  expect_equal(info$testing$effective_N, 1L)
  expect_equal(info$testing$effective_n, 1L)
  expect_equal(info$testing$tested_count, 1L)
})

test_that("background diagnostics reconcile invalid, repeated, and unmapped IDs", {
  annot <- data.frame(GeneID = c("A", "B"), Term = c("T", "T"))
  result <- suppressWarnings(suppressMessages(richR::enrich(
    c("A", "A", "B", "X", NA, ""), annot,
    universe = c("A", "A", "Y", NA, " "),
    pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  )))
  info <- richR::enrichmentInfo(result)
  expect_equal(unname(unlist(info$background[c("original_count", "valid_count", "unique_count",
                                               "invalid_count", "duplicate_count", "unmapped_count",
                                               "effective_count")])),
               c(5L, 3L, 2L, 2L, 1L, 1L, 1L))
  expect_equal(info$background$unmapped_ids, "Y")
  expect_equal(info$background$effective_ids, "A")
})

test_that("enrichmentInfo handles legacy objects without inventing provenance", {
  annot <- new("Annot", species = "old", anntype = "custom",
               keytype = "SYMBOL", annot = data.frame(GeneID = "A", Term = "T"))
  legacy <- legacy_rich_result()
  gsea <- legacy_gsea_result()
  expected <- list(schema_version = 1L, available = FALSE,
                   reason = "metadata not recorded", diagnostics = NULL,
                   provenance = NULL)
  expect_equal(richR::enrichmentInfo(annot), expected)
  expect_equal(richR::enrichmentInfo(legacy), expected)
  expect_equal(richR::enrichmentInfo(gsea), expected)
  legacy_roundtrips <- list()
  for (i in seq_along(list(annot, legacy, gsea))) {
    object <- list(annot, legacy, gsea)[[i]]
    restored_path <- tempfile(fileext = ".rds")
    saveRDS(object, restored_path)
    restored <- readRDS(restored_path)
    legacy_roundtrips[[i]] <- restored
    expect_true(methods::validObject(restored))
    expect_equal(richR::enrichmentInfo(restored), expected)
  }
  expect_s3_class(richR::result(legacy_roundtrips[[2]]), "data.frame")
  expect_s3_class(richR::result(legacy_roundtrips[[3]]), "data.frame")
  expect_error(richR::enrichmentInfo(data.frame()), "Annot, richResult, or GSEAResult")
})

test_that("annotation constructors record custom source and normalization history", {
  annot <- suppressWarnings(suppressMessages(richR::buildAnnotFromList(
    list(T = c("A", "A", NA_character_, " ", "B"))
  )))
  info <- richR::enrichmentInfo(annot)
  expect_true(info$available)
  expect_equal(info$analysis, "annotation")
  expect_equal(unname(unlist(info$annotation[c("original_rows", "valid_rows", "invalid_rows",
                                               "duplicate_rows", "unique_rows")])),
               c(5L, 3L, 2L, 1L, 2L))
  expect_equal(info$provenance$annotation$source_type, "custom")
  expect_equal(info$provenance$annotation$constructor, "buildAnnotFromList")
  expect_equal(info$provenance$annotation$constructor_normalization$duplicate_rows, 1L)
  annot_path <- tempfile(fileext = ".rds")
  saveRDS(annot, annot_path)
  annot_restored <- readRDS(annot_path)
  expect_true(methods::validObject(annot_restored))
  expect_equal(richR::enrichmentInfo(annot_restored), info)

  gmt_file <- tempfile(fileext = ".gmt")
  writeLines("T\tdescription\tA\tB", gmt_file)
  gmt <- suppressMessages(richR::readGMT(gmt_file))
  gmt_info <- richR::enrichmentInfo(gmt)
  expect_equal(gmt_info$provenance$annotation$source_type, "custom")
  expect_equal(gmt_info$provenance$annotation$source_name, basename(gmt_file))
  expect_true(is.na(gmt_info$provenance$annotation$source_version))
})

test_that("package-backed constructors retain observed provider provenance", {
  go_rows <- data.frame(GeneID = c("A", NA_character_), GOALL = c("GO:1", "GO:2"),
                        ONTOLOGYALL = c("BP", "BP"), Annot = c("term", "term"))
  go <- expect_warning(testthat::with_mocked_bindings(
    .makeGOdata = function(...) go_rows,
    .getdbname = function(...) "stats",
    richR::buildAnnot(species = "human", keytype = "SYMBOL", anntype = "GO"),
    .package = "richR"
  ), "excluded 1")
  go_info <- richR::enrichmentInfo(go)
  expect_equal(go_info$provenance$annotation$provider_package, "stats")
  expect_equal(go_info$provenance$annotation$provider_version,
               as.character(utils::packageVersion("stats")))
  expect_equal(go_info$provenance$annotation$constructor_normalization$invalid_rows, 1L)

  msig_rows <- data.frame(gs_cat = "H", gs_subcat = NA_character_,
                          gene_symbol = c("A", NA_character_), gs_exact_source = c("T", "T"),
                          gs_name = c("HALLMARK_T", "HALLMARK_T"),
                          db_version = "2025.1")
  msig <- expect_warning(testthat::with_mocked_bindings(
    .getmsig = function(...) "Homo sapiens",
    msigdbr = function(...) msig_rows,
    richR::buildMSIGDB(anntype = "HALLMARK"),
    .package = "richR"
  ), "excluded 1")
  msig_info <- richR::enrichmentInfo(msig)
  expect_equal(msig_info$provenance$annotation$source_version, "2025.1")
  expect_equal(msig_info$provenance$annotation$provider_package, "msigdbr")
  expect_equal(msig_info$provenance$annotation$provider_version,
               as.character(utils::packageVersion("msigdbr")))

  own <- expect_warning(testthat::with_mocked_bindings(
    keys = function(x, keytype) "A",
    getann = function(...) data.frame(annotation = c("term", "other"),
                                      row.names = c("GO:1", "GO:2")),
    .package = "richR",
    testthat::with_mocked_bindings(
      select = function(x, keys, keytype, columns) {
        data.frame(GeneID = c("A", NA_character_), GOALL = c("GO:1", "GO:2"),
                   ONTOLOGYALL = c("BP", "BP"))
      },
      richR::buildOwn(dbname = "org.Hs.eg.db", anntype = "GO", species = "mouse", keytype = "SYMBOL"),
      .package = "AnnotationDbi"
    )
  ), "excluded 1")
  own_info <- richR::enrichmentInfo(own)
  expect_equal(own_info$provenance$annotation$source_name, "org.Hs.eg.db")
  expect_equal(own_info$provenance$annotation$provider_version,
               as.character(utils::packageVersion("org.Hs.eg.db")))
  expect_equal(own_info$provenance$annotation$constructor_normalization$invalid_rows, 1L)
})

test_that("live KEGG annotation records KEGGREST as its mapping provider", {
  kegg <- testthat::with_mocked_bindings(
    .getdbname = function(...) "org.Hs.eg.db",
    .getspeices = function(...) "hsa",
    keggLink = function(...) c("hsa:1" = "path:hsa00010", "hsa:2" = "path:hsa00020"),
    getann = function(...) data.frame(annotation = c("glycolysis", "citrate cycle"),
                                      row.names = c("00010", "00020")),
    richR::buildAnnot(species = "human", keytype = "ENTREZID", anntype = "KEGG",
                      builtin = FALSE),
    .package = "richR"
  )
  info <- richR::enrichmentInfo(kegg)
  expect_equal(info$provenance$annotation$source_type, "KEGG")
  expect_equal(info$provenance$annotation$source_name, "KEGG")
  expect_equal(info$provenance$annotation$provider_package, "KEGGREST")
  expect_equal(info$provenance$annotation$provider_version,
               as.character(utils::packageVersion("KEGGREST")))
  expect_true(is.na(info$provenance$annotation$source_version))
})

test_that("KEGG module constructor records its actual membership provider", {
  module <- testthat::with_mocked_bindings(
    .getdbname = function(...) "org.Hs.eg.db",
    .getspeices = function(...) "hsa",
    keggLink = function(...) c("hsa:1" = "md:hsa:M00001",
                               "hsa:2" = "md:hsa:M00002"),
    richR::buildAnnot(species = "human", keytype = "ENTREZID", anntype = "KEGGM"),
    .package = "richR"
  )
  info <- richR::enrichmentInfo(module)
  expect_equal(info$provenance$annotation$source_type, "KEGG")
  expect_equal(info$provenance$annotation$source_name, "KEGG")
  expect_equal(info$provenance$annotation$provider_package, "KEGGREST")
  expect_equal(info$provenance$annotation$provider_version,
               as.character(utils::packageVersion("KEGGREST")))
  expect_true(is.na(info$provenance$annotation$source_version))
})

test_that("GO constructor normalizes duplicate rows from the actual helper path", {
  go <- testthat::with_mocked_bindings(
    select = function(x, keys, keytype, columns) {
      data.frame(GeneID = c("1", "1", "2"), GOALL = "GO:1",
                 ONTOLOGYALL = "BP")
    },
    .package = "AnnotationDbi",
    testthat::with_mocked_bindings(
      .getdbname = function(...) "org.Hs.eg.db",
      keys = function(x, keytype) c("1", "2"),
      getann = function(...) data.frame(annotation = "term", row.names = "GO:1"),
      richR::buildAnnot(species = "human", keytype = "ENTREZID", anntype = "GO"),
      .package = "richR"
    )
  )
  info <- richR::enrichmentInfo(go)
  expect_equal(unname(unlist(info$provenance$annotation$constructor_normalization[
    c("original_rows", "valid_rows", "duplicate_rows", "unique_rows")])),
    c(3L, 3L, 1L, 2L))
})

test_that("analysis metadata survives serialization and object-preserving transforms", {
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("T1", "T1", "T2", "T2"))
  result <- suppressMessages(richR::enrich(
    c("A", "B"), annot, pvalue = 1, minSize = 1,
    minGSSize = 1, keepRich = FALSE
  ))
  original_info <- richR::enrichmentInfo(result)
  expect_true(is.list(attr(result, "richR_metadata")))
  rds <- tempfile(fileext = ".rds")
  saveRDS(result, rds)
  restored <- readRDS(rds)
  expect_s4_class(restored, "richResult")
  expect_equal(richR::enrichmentInfo(restored), original_info)
  filtered <- dplyr::filter(result, Pvalue <= 1)
  arranged <- dplyr::arrange(result, dplyr::desc(Pvalue))
  expect_equal(richR::enrichmentInfo(filtered), original_info)
  expect_equal(richR::enrichmentInfo(arranged), original_info)
  extracted <- result[1, , drop = FALSE]
  expect_true(is.data.frame(extracted))
  expect_null(attr(extracted, "richR_metadata"))
})

test_that("data-frame detail retains caller columns while metadata tracks row names", {
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("T", "T", "U", "U"))
  query <- data.frame(logFC = c(2, 3), row.names = c("A", "B"))
  result <- suppressMessages(richR::enrich(
    query, annot, pvalue = 1, minSize = 1, minGSSize = 1, keepRich = FALSE
  ))
  expect_true("logFC" %in% names(result@detail))
  expect_equal(result@detail$logFC, c(2, 3))
  expect_equal(richR::enrichmentInfo(result)$input$ids, c("A", "B"))
})

test_that("public ORA empty reports retain complete testing diagnostics", {
  annot <- data.frame(GeneID = "A", Term = "T")
  result <- suppressMessages(richR::enrich(
    "A", annot, pvalue = 1, minSize = 2, minGSSize = 1,
    keepRich = FALSE
  ))
  expect_s4_class(result, "richResult")
  expect_equal(nrow(result@result), 0L)
  expect_equal(nrow(result@detail), 0L)
  info <- richR::enrichmentInfo(result)
  expect_equal(info$testing$effective_N, 1L)
  expect_equal(info$testing$effective_n, 1L)
  expect_equal(info$testing$tested_count, 1L)
})

test_that("GSEA metadata describes the tested backend family and ranked input", {
  annot <- data.frame(GeneID = c("A", "B", "C", "D"),
                      Term = c("T1", "T1", "T2", "T2"))
  previous <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(previous), add = TRUE)
  result <- richR::richGSEA(
    c(A = 2, B = 1, C = -1, D = -2), annot, minSize = 1,
    maxSize = 5, pvalue = 1, table = FALSE
  )
  info <- richR::enrichmentInfo(result)
  expect_equal(info$analysis, "GSEA")
  expect_equal(info$input$ids, c("A", "B", "C", "D"))
  expect_equal(info$background$mode, "ranked_input")
  expect_equal(info$background$effective_ids, c("A", "B", "C", "D"))
  expect_equal(info$testing$tested_count, 2L)
  expect_equal(info$testing$family_policy, "fgsea_returned_before_filtering")
  expect_equal(info$provenance$fgsea_version, as.character(utils::packageVersion("fgsea")))
  expect_equal(info$annotation$effective_rows, info$annotation$eligible_rows)
  gsea_restored <- tempfile(fileext = ".rds")
  saveRDS(result, gsea_restored)
  expect_true(methods::validObject(readRDS(gsea_restored)))
  expect_equal(richR::enrichmentInfo(readRDS(gsea_restored)), info)
})
