gsea_test_annotation <- function() {
  data.frame(GeneID = c("A", "B", "C", "D"),
             Term = c("T1", "T1", "T2", "T2"))
}

gsea_recording_backend <- function(...) {
  args <- list(...)
  if (!exists(".gsea_backend_args", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".gsea_backend_args", list(), envir = .GlobalEnv)
  }
  assign(".gsea_backend_args", c(get(".gsea_backend_args", envir = .GlobalEnv,
                                      inherits = FALSE), list(args)),
         envir = .GlobalEnv)
  data.frame(pathway = "T1", pval = 0.01, NES = 1, size = 2,
             leadingEdge = I(list(c("A", "B"))))
}

test_that("GSEA controls are explicitly omitted or forwarded through both methods", {
  old_args <- if (exists(".gsea_backend_args", envir = .GlobalEnv, inherits = FALSE)) {
    get(".gsea_backend_args", envir = .GlobalEnv)
  } else NULL
  on.exit(if (is.null(old_args)) {
    if (exists(".gsea_backend_args", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".gsea_backend_args", envir = .GlobalEnv)
    }
  } else assign(".gsea_backend_args", old_args, envir = .GlobalEnv), add = TRUE)
  assign(".gsea_backend_args", list(), envir = .GlobalEnv)
  annot <- gsea_test_annotation()
  stats <- c(A = 2, B = 1, C = -1, D = -2)

  result <- testthat::with_mocked_bindings(
    fgseaMultilevel = gsea_recording_backend,
    richR::richGSEA(stats, annot, minSize = 1),
    .package = "richR"
  )
  expect_s4_class(result, "GSEAResult")
  args <- get(".gsea_backend_args", envir = .GlobalEnv)[[1]]
  expect_setequal(names(args), c("pathways", "stats", "minSize", "maxSize"))

  testthat::with_mocked_bindings(
    fgseaMultilevel = gsea_recording_backend,
    richR::richGSEA(stats, annot, minSize = 1, maxSize = Inf),
    .package = "richR"
  )
  args <- tail(get(".gsea_backend_args", envir = .GlobalEnv), 1)[[1]]
  expect_identical(args$maxSize, Inf)

  positional <- testthat::with_mocked_bindings(
    fgseaMultilevel = gsea_recording_backend,
    richR::richGSEA(stats, annot, "SYMBOL", 1, NULL, FALSE, 1, "", Inf,
                    "BH", NULL, FALSE, "|"),
    .package = "richR"
  )
  expect_identical(positional@sep, "|")
  expect_identical(positional@pvalueCutoff, 1)

  controls <- list(eps = 0, scoreType = "pos", nproc = 1,
                   BPPARAM = BiocParallel::SerialParam())
  for (control_name in names(controls)) {
    call_args <- c(list(stats, annot, minSize = 1), controls[control_name])
    testthat::with_mocked_bindings(
      fgseaMultilevel = gsea_recording_backend,
      do.call(richR::richGSEA, call_args),
      .package = "richR"
    )
    args <- tail(get(".gsea_backend_args", envir = .GlobalEnv), 1)[[1]]
    expect_identical(args[[control_name]], controls[[control_name]])
    expect_length(intersect(names(args), setdiff(names(controls), control_name)), 0L)
  }

  annot_object <- new("Annot", species = "custom", anntype = "custom",
                      keytype = "SYMBOL", annot = annot)
  testthat::with_mocked_bindings(
    fgseaMultilevel = gsea_recording_backend,
    richR::richGSEA(stats, annot_object, minSize = 1, eps = 0, scoreType = "pos"),
    .package = "richR"
  )
  args <- tail(get(".gsea_backend_args", envir = .GlobalEnv), 1)[[1]]
  expect_identical(args$eps, 0)
  expect_identical(args$scoreType, "pos")
  expect_false("nproc" %in% names(args))
  expect_false("BPPARAM" %in% names(args))
})

test_that("GSEA rejects malformed controls and ranked vectors before backend execution", {
  .gsea_backend_called <- FALSE
  invalid_controls <- list(
    list(eps = numeric()), list(eps = c(0, 1)), list(eps = NA_real_),
    list(eps = Inf), list(eps = -0.1), list(scoreType = "bad"),
    list(scoreType = c("pos", "neg")), list(scoreType = NA_character_),
    list(scoreType = factor("pos")),
    list(nproc = -1), list(nproc = 1.5), list(nproc = Inf),
    list(nproc = NA_real_), list(BPPARAM = "serial"),
    list(nproc = 1, BPPARAM = BiocParallel::SerialParam()),
    list(seed = -1), list(seed = 1.5), list(seed = Inf), list(seed = NA_real_),
    list(unknown_backend_option = TRUE)
  )
  for (control in invalid_controls) {
    expect_error(testthat::with_mocked_bindings(
      fgseaMultilevel = function(...) {
        .gsea_backend_called <<- TRUE
        gsea_recording_backend(...)
      },
      do.call(richR::richGSEA,
              c(list(c(A = 1, B = -1), gsea_test_annotation(), minSize = 1), control)),
      .package = "richR"
    ))
  }
  expect_false(.gsea_backend_called)

  invalid_inputs <- list(c(1, 2), c(A = 1, A = -1), c("bad", "score"),
                         c(A = NA_real_, B = 1), c(A = Inf, B = 1),
                         c(A = 1, B = -1)[FALSE])
  for (input in invalid_inputs) {
    expect_error(richR:::richGSEA_internal(input, gsea_test_annotation(), minSize = 1),
                 "named numeric vector|unique|finite|non-empty")
  }
  expect_error(richR::richGSEA(c(A = 1, B = -1), gsea_test_annotation(),
                               minSize = 1, unsupported = TRUE),
               "unused argument|unused.*unsupported")
})

test_that("parGSEA forwards controls and preserves both return modes", {
  old_args <- if (exists(".gsea_backend_args", envir = .GlobalEnv, inherits = FALSE)) {
    get(".gsea_backend_args", envir = .GlobalEnv)
  } else NULL
  on.exit(if (is.null(old_args)) {
    if (exists(".gsea_backend_args", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".gsea_backend_args", envir = .GlobalEnv)
    }
  } else assign(".gsea_backend_args", old_args, envir = .GlobalEnv), add = TRUE)
  assign(".gsea_backend_args", list(), envir = .GlobalEnv)
  annot <- gsea_test_annotation()
  data <- data.frame(Gene = c("A", "B", "C", "D"),
                     log2FoldChange = c(2, 1, -1, -2))
  rownames(data) <- data$Gene

  result <- testthat::with_mocked_bindings(
    fgseaMultilevel = gsea_recording_backend,
    richR::parGSEA(data, annot, log2FC = "log2FoldChange", eps = 0,
                   simplify = FALSE),
    .package = "richR"
  )
  expect_s4_class(result, "GSEAResult")
  expect_equal(result@input, c(D = -2, C = -1, B = 1, A = 2))
  expect_true(is.list(attr(result, "richR_metadata")))
  expect_identical(tail(get(".gsea_backend_args", envir = .GlobalEnv), 1)[[1]]$eps, 0)

  result <- testthat::with_mocked_bindings(
    fgseaMultilevel = gsea_recording_backend,
    richR::parGSEA(data, annot, gene.col = "Gene", simplify = TRUE,
                   scoreType = "neg"),
    .package = "richR"
  )
  expect_s3_class(result, "data.frame")
  expect_identical(tail(get(".gsea_backend_args", envir = .GlobalEnv), 1)[[1]]$scoreType, "neg")
})

test_that("seeded serial GSEA restores RNG and reproduces results", {
  old_param <- BiocParallel::bpparam()
  BiocParallel::register(BiocParallel::SerialParam())
  on.exit(BiocParallel::register(old_param), add = TRUE)
  annot <- data.frame(GeneID = rep(paste0("g", 1:20), 2),
                      Term = rep(c("left", "right"), each = 20))
  stats <- setNames(seq(-2, 2, length.out = 20), paste0("g", 1:20))
  set.seed(928)
  old_kind <- RNGkind()
  old_seed <- .Random.seed
  old_kind <- RNGkind()
  first <- suppressMessages(richR::richGSEA(
    stats, annot, minSize = 5, maxSize = 30, pvalue = 1, seed = 77
  ))
  expect_identical(RNGkind(), old_kind)
  expect_identical(.Random.seed, old_seed)
  second <- suppressMessages(richR::richGSEA(
    stats, annot, minSize = 5, maxSize = 30, pvalue = 1, seed = 77
  ))
  expect_equal(first@result, second@result)
  expect_equal(richR::enrichmentInfo(first)$provenance$rng$seed, 77L)
  expect_false(".Random.seed" %in% names(richR::enrichmentInfo(first)$provenance$rng))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) .Random.seed else NULL
  rm(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv), add = TRUE)
  suppressMessages(richR::richGSEA(
    stats, annot, minSize = 5, maxSize = 30, pvalue = 1, seed = 77
  ))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), old_kind)
})

test_that("seeded GSEA restores RNG after backend errors and seed NULL advances it", {
  set.seed(54)
  old_kind <- RNGkind()
  old_seed <- .Random.seed
  expect_error(testthat::with_mocked_bindings(
    fgseaMultilevel = function(...) stop("backend failed"),
    richR::richGSEA(c(A = 1, B = -1), gsea_test_annotation(), minSize = 1, seed = 7),
    .package = "richR"
  ), "backend failed")
  expect_identical(RNGkind(), old_kind)
  expect_identical(.Random.seed, old_seed)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved_seed <- if (had_seed) .Random.seed else NULL
  rm(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv), add = TRUE)
  expect_error(testthat::with_mocked_bindings(
    fgseaMultilevel = function(...) stop("backend failed"),
    richR::richGSEA(c(A = 1, B = -1), gsea_test_annotation(), minSize = 1, seed = 7),
    .package = "richR"
  ), "backend failed")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_identical(RNGkind(), old_kind)

  set.seed(66)
  old_seed <- .Random.seed
  testthat::with_mocked_bindings(
    fgseaMultilevel = function(...) {
      runif(1)
      gsea_recording_backend(...)
    },
    richR::richGSEA(c(A = 1, B = -1), gsea_test_annotation(), minSize = 1),
    .package = "richR"
  )
  expect_false(identical(.Random.seed, old_seed))
})
