.run_consumer_rscript <- function(code) {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  lib_paths <- paste(capture.output(dput(.libPaths())), collapse = "")
  writeLines(c(
    paste0(".libPaths(", lib_paths, ")"),
    "setwd(tempdir())",
    code
  ), script)
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  list(status = if (is.null(status)) 0L else status, output = output)
}

.expect_consumer_process_success <- function(run) {
  cat(run$output, sep = "\n")
  expect_identical(run$status, 0L, info = paste(run$output, collapse = "\n"))
}

.cluster_fixture <- function() {
  data.frame(
    Term = LETTERS[1:5], Annot = letters[1:5],
    GeneID = c(
      paste0("g", 1:5, collapse = ","),
      paste0("g", 2:6, collapse = ","),
      paste0("g", 11:15, collapse = ","),
      paste0("g", 12:16, collapse = ","),
      paste0("g", 17:40, collapse = ",")
    ),
    Pvalue = c(1e-8, 1e-7, 1e-6, 1e-5, .5),
    Padj = c(1e-8, 1e-7, 1e-6, 1e-5, .5)
  )
}

.as_rich_result <- function(result) {
  loadNamespace("richR")
  methods::new(
    "richResult", result = result, detail = data.frame(),
    pvalueCutoff = .1, pAdjustMethod = "BH", padjCutoff = .1,
    genenumber = length(unique(unlist(strsplit(result$GeneID, ",", fixed=TRUE)))),
    organism = "", ontology = "",
    gene = unique(unlist(strsplit(result$GeneID, ",", fixed=TRUE))),
    keytype = "", sep = ","
  )
}

test_that("namespace-only KEGG enrichment resolves bundled data and keeps metadata", {
  run <- .run_consumer_rscript(paste(
    "stopifnot(requireNamespace('richR', quietly=TRUE))",
    "stopifnot(!('package:richR' %in% search()))",
    "cat('PACKAGE_VERSION:', as.character(utils::packageVersion('richR')), '\\n')",
    "package_path <- getNamespaceInfo(asNamespace('richR'), 'path'); cat('PACKAGE_PATH:', package_path, '\\n')",
    "cat('WORKING_DIRECTORY:', getwd(), '\\n')",
    "stopifnot(!startsWith(normalizePath(getwd()), normalizePath(package_path)))",
    "ref <- new.env(parent=emptyenv()); utils::data('path', package='richR', envir=ref)",
    "ids <- unique(as.character(ref$path$ko)); ids <- ids[!is.na(ids) & nzchar(ids)]",
    "annot <- data.frame(GeneID=c('A','B','C','D'), Term=rep(ids[1:2], each=2))",
    "warnings <- character()",
    "res <- withCallingHandlers(richR::richKEGG(c('A','B'), annot, pvalue=1, minSize=1, minGSSize=1, keepRich=FALSE), warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') })",
    "stopifnot(!any(grepl('data set.*not found', warnings, ignore.case=TRUE)))",
    "stopifnot(nrow(res@result)>0L, identical(as.character(res@result$ko), '00010'))",
    "stopifnot(all(c('ko','Level1','Level2','Level3') %in% names(res@result)))",
    "expected <- ref$path[match(res@result$ko, ref$path$ko), c('ko','Level3','Level2','Level1'), drop=FALSE]",
    "stopifnot(identical(as.character(res@result$ko), as.character(expected$ko)))",
    "stopifnot(identical(as.character(res@result$Level3), as.character(expected$Level3)))",
    "stopifnot(identical(as.character(res@result$Level2), as.character(expected$Level2)))",
    "stopifnot(identical(as.character(res@result$Level1), as.character(expected$Level1)))",
    "empty <- richR::richKEGG(c('A','B'), annot, pvalue=1e-10, minSize=1, minGSSize=1, keepRich=FALSE)",
    "stopifnot(is.data.frame(empty@result), nrow(empty@result)==0L, identical(names(empty@result), names(res@result)))",
    "base::library(richR)",
    "attached <- richR::richKEGG(c('A','B'), annot, pvalue=1, minSize=1, minGSSize=1, keepRich=FALSE)",
    "base::detach('package:richR', unload=FALSE)",
    "stopifnot(isTRUE(all.equal(res@result, attached@result, check.attributes=FALSE)))",
    "stopifnot(!('package:richR' %in% search()))",
    "cat('NAMESPACE_ONLY_KEGG_PASS\\n')",
    sep = "\n"
  ))
  .expect_consumer_process_success(run)
  expect_true(any(grepl("NAMESPACE_ONLY_KEGG_PASS", run$output)))
  expect_true(any(grepl("WORKING_DIRECTORY:", run$output)))
})

test_that("bundled data access preserves global state and works in a worker process", {
  run <- .run_consumer_rscript(paste(
    "stopifnot(requireNamespace('richR', quietly=TRUE), !('package:richR' %in% search()))",
    "cat('PACKAGE_VERSION:', as.character(utils::packageVersion('richR')), '\\n')",
    "package_path <- getNamespaceInfo(asNamespace('richR'), 'path'); cat('PACKAGE_PATH:', package_path, '\\n')",
    "cat('WORKING_DIRECTORY:', getwd(), '\\n')",
    "stopifnot(!startsWith(normalizePath(getwd()), normalizePath(package_path)))",
    "ref <- new.env(parent=emptyenv()); utils::data('kegg', package='richR', envir=ref); expected_kegg <- ref$kegg.db",
    "ref <- new.env(parent=emptyenv()); utils::data('module', package='richR', envir=ref); expected_module <- ref$module",
    "marker <- data.frame(sentinel=TRUE); path <- marker; kegg.db <- marker; module <- marker",
    "stopifnot(identical(getFromNamespace('.get_kg_dat','richR')(), expected_kegg))",
    "stopifnot(identical(getFromNamespace('.get_kgm.data','richR')(), expected_module))",
    "stopifnot(identical(path,marker), identical(kegg.db,marker), identical(module,marker))",
    "rm(path,kegg.db,module)",
    "path_ids <- c('00010','00020')",
    "path_ref <- new.env(parent=emptyenv()); utils::data('path', package='richR', envir=path_ref)",
    "level_annot <- data.frame(GeneID=c('A','B','C','D'), PATH=rep(path_ids, each=2))",
    "level_warnings <- character()",
    "level_result <- withCallingHandlers(richR::richLevel(c('A','B'), level_annot, level='Level3', pvalue=1, minSize=1, minGSSize=1, keepRich=FALSE), warning=function(w) { level_warnings <<- c(level_warnings, conditionMessage(w)); invokeRestart('muffleWarning') })",
    "stopifnot(!any(grepl('data set.*not found', level_warnings, ignore.case=TRUE)), methods::is(level_result,'richResult'), nrow(level_result@result)>0L)",
    "expected_level <- path_ref$path$Level3[match(level_result@result$Term, path_ref$path$Level3)]",
    "stopifnot(identical(as.character(level_result@result$Term), as.character(expected_level)))",
    "scores <- setNames(c(-10:-1,1:10), paste0('g',1:20))",
    "gsea_annot <- data.frame(GeneID=paste0('g',c(1:6,7:12)), Term=rep(path_ref$path$Level3[1:2], each=6))",
    "gsea_warnings <- character()",
    "set.seed(42)",
    "gsea_result <- withCallingHandlers(richR::richGSEA(scores, object=gsea_annot, KEGG=TRUE, minSize=1, maxSize=20, pvalue=1, padj=1, organism='hsa'), warning=function(w) { gsea_warnings <<- c(gsea_warnings, conditionMessage(w)); invokeRestart('muffleWarning') })",
    "stopifnot(!any(grepl('data set.*not found', gsea_warnings, ignore.case=TRUE)), methods::is(gsea_result,'GSEAResult'))",
    "stopifnot(all(c('Level1','Level2','Level3') %in% names(gsea_result@result)))",
    "stopifnot(nrow(gsea_result@result)==2L)",
    "stopifnot(setequal(as.character(gsea_result@result$pathway), as.character(path_ref$path$Level3[1:2])))",
    "expected_gsea <- path_ref$path[match(gsea_result@result$pathway, path_ref$path$Level3), c('Level3','Level2','Level1'), drop=FALSE]",
    "stopifnot(identical(as.character(gsea_result@result$Level3), as.character(expected_gsea$Level3)))",
    "stopifnot(identical(as.character(gsea_result@result$Level2), as.character(expected_gsea$Level2)))",
    "stopifnot(identical(as.character(gsea_result@result$Level1), as.character(expected_gsea$Level1)))",
    "stopifnot(identical(getFromNamespace('.get_kg_dat','richR')(), expected_kegg))",
    "stopifnot(identical(getFromNamespace('.get_kgm.data','richR')(), expected_module))",
    "stopifnot(!exists('path',.GlobalEnv,inherits=FALSE), !exists('kegg.db',.GlobalEnv,inherits=FALSE), !exists('module',.GlobalEnv,inherits=FALSE))",
    "stopifnot(!('package:richR' %in% search()))",
    "cat('NAMESPACE_ONLY_GLOBAL_BINDINGS_PASS\\n')",
    sep = "\n"
  ))
  .expect_consumer_process_success(run)
  expect_true(any(grepl("NAMESPACE_ONLY_GLOBAL_BINDINGS_PASS", run$output)))
  expect_true(any(grepl("WORKING_DIRECTORY:", run$output)))
})

test_that("balanced cluster memberships and scores align for both public representations", {
  fixture <- .cluster_fixture()
  for (input in list(fixture, .as_rich_result(fixture))) {
    for (use_term in c(FALSE, TRUE)) {
      result <- tryCatch(richR::richCluster(
        input, gene=NULL, useTerm=use_term, cutoff=.1, overlap=0, minSize=2
      ), error=function(e) e)
      expect_false(inherits(result, "error"), info=if (inherits(result, "error")) conditionMessage(result))
      if (inherits(result, "error")) next
      expected <- if (use_term) c("A", "C") else c("a", "c")
      expected_members <- if (use_term) list(c("A", "B"), c("C", "D")) else
        list(c("a", "b"), c("c", "d"))
      expect_identical(rownames(result), expected)
      expect_identical(as.integer(result$filteredClusterSize), c(2L, 2L))
      expect_equal(result$EnrichmentScore, c(7.5, 5.5), tolerance=1e-12)
      expect_true(all(vapply(seq_len(2), function(i) {
        setequal(strsplit(result$Cluster[i], ",", fixed=TRUE)[[1]], expected_members[[i]])
      }, logical(1))))
      filtered <- richR::richCluster(
        input, gene=NULL, useTerm=use_term, cutoff=.1, overlap=0,
        minSize=2, escore=6
      )
      expect_identical(rownames(filtered), expected[1])
      expect_identical(as.integer(filtered$filteredClusterSize), 2L)
      expect_equal(filtered$EnrichmentScore, 7.5, tolerance=1e-12)
    }
  }
})

test_that("regular three-term cliques survive initial neighbor grouping", {
  clique <- data.frame(
    Term=c("A", "B", "C", "E"), Annot=c("a", "b", "c", "e"),
    GeneID=c("g1,g2,g3,g4,g5", "g1,g2,g3,g4,g5", "g1,g2,g3,g4,g5",
             paste0("g", 17:40, collapse=",")),
    Pvalue=c(1e-8, 1e-7, 1e-6, .5), Padj=c(1e-8, 1e-7, 1e-6, .5)
  )
  for (input in list(clique, .as_rich_result(clique))) {
    for (use_term in c(FALSE, TRUE)) {
      result <- tryCatch(richR::richCluster(
        input, gene=NULL, useTerm=use_term, cutoff=.1, overlap=0, minSize=2
      ), error=function(e) e)
      expect_false(inherits(result, "error"), info=if (inherits(result, "error")) conditionMessage(result))
      if (inherits(result, "error")) next
      expect_identical(rownames(result), if (use_term) "A" else "a")
      expect_identical(as.integer(result$filteredClusterSize), 3L)
      expect_equal(result$EnrichmentScore, 7, tolerance=1e-12)
      expect_setequal(strsplit(result$Cluster, ",", fixed=TRUE)[[1]],
                      if (use_term) c("A", "B", "C") else c("a", "b", "c"))
    }
  }
})

test_that("unequal clusters preserve membership, scores, and strict filtering", {
  unequal <- data.frame(
    Term=LETTERS[1:6], Annot=letters[1:6],
    GeneID=c(paste0("g", 1:5, collapse=","), paste0("g", 2:6, collapse=","),
             paste0("g", 3:7, collapse=","), paste0("g", 11:15, collapse=","),
             paste0("g", 12:16, collapse=","), paste0("g", 17:40, collapse=",")),
    Pvalue=c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, .5),
    Padj=c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, .5)
  )
  result <- richR::richCluster(unequal, gene=NULL, cutoff=.1, overlap=0, minSize=2)
  expect_identical(rownames(result), c("a", "d"))
  expect_identical(as.integer(result$filteredClusterSize), c(3L, 2L))
  expect_equal(result$EnrichmentScore, c(7, 4.5), tolerance=1e-12)
  expect_setequal(strsplit(result$Cluster[1], ",", fixed=TRUE)[[1]], letters[1:3])
  expect_setequal(strsplit(result$Cluster[2], ",", fixed=TRUE)[[1]], letters[4:5])

  filtered <- richR::richCluster(
    unequal, gene=NULL, cutoff=.1, overlap=0, minSize=2, escore=6
  )
  expect_identical(rownames(filtered), "a")
  expect_identical(as.integer(filtered$filteredClusterSize), 3L)
})

test_that("cluster errors remain clear when no pairs or scores survive", {
  fixture <- .cluster_fixture()
  expect_error(
    richR::richCluster(fixture, gene=NULL, cutoff=1, overlap=0, minSize=2),
    "kappa cutoff too high"
  )
  expect_error(
    richR::richCluster(fixture, gene=NULL, cutoff=.1, overlap=0, minSize=2, escore=8),
    "Enrichment score too high, No significant cluster"
  )
})
