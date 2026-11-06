# Result objects without a recorded tested table (SPEC-RICHR-004 REQ-005,
# AC-004): objects serialized by richR before the tested slot existed
# (fixtures/legacy_*_fca5276.rds, made by fixtures/make-legacy-fixtures.R with
# the fca5276 build), objects built by as(), and objects whose slot and
# metadata are gone. Catches a tested table fabricated from @result, metadata
# silently upgraded, and any code path that reads x@tested without checking
# for it or otherwise changes what pre-0.2.0 objects return.

unavailable_record <- list(schema_version = 1L, available = FALSE,
                           reason = "metadata not recorded",
                           diagnostics = NULL, provenance = NULL)

legacy_fixture <- function(name) readRDS(test_path("fixtures", name))

test_that("a pre-0.2.0 richResult has no tested table and returns what fca5276 returned", {
  x <- legacy_fixture("legacy_richResult_fca5276.rds")
  ref <- legacy_fixture("legacy_reference_fca5276.rds")$richResult
  expect_false(methods::.hasSlot(x, "tested"))
  expect_gt(nrow(x@result), 0)

  expect_warning(tt <- testedTerms(x), "tested table unavailable: .*no tested slot")
  expect_null(tt)

  info <- enrichmentInfo(x)
  expect_identical(info, attr(x, "richR_metadata"))
  expect_identical(info$schema_version, 1L)
  expect_null(info$testing$zero_overlap_count)

  expect_identical(result(x), ref$result)
  expect_identical(compareResult(list(a = x)), ref$compareResult)
  view <- dplyr::filter(x, Pvalue <= 0.5)
  expect_identical(view@result, ref$filter)
  expect_false(methods::.hasSlot(view, "tested"))
  # libm log10 differs by one bit between Linux and Windows. Keep the stored
  # analysis frame exact, and compute only the derived plotting column locally.
  ref$richBar$neg_log_p <- -log10(as.numeric(ref$richBar$Padj))
  expect_identical(richBar(x)$data, ref$richBar)
})

test_that("a pre-0.2.0 GSEAResult has no tested table and returns what fca5276 returned", {
  x <- legacy_fixture("legacy_GSEAResult_fca5276.rds")
  ref <- legacy_fixture("legacy_reference_fca5276.rds")$GSEAResult
  expect_false(methods::.hasSlot(x, "tested"))
  expect_gt(nrow(x@result), 0)

  expect_warning(tt <- testedTerms(x), "tested table unavailable: .*no tested slot")
  expect_null(tt)

  info <- enrichmentInfo(x)
  expect_identical(info, attr(x, "richR_metadata"))
  expect_identical(info$schema_version, 1L)

  expect_identical(result(x), ref$result)
  view <- dplyr::filter(x, pval <= 0.5)
  expect_identical(view@result, ref$filter)
  expect_false(methods::.hasSlot(view, "tested"))
})

test_that("an as()-built richResult has no tested table and reads its own @result", {
  l <- as(.fora_reference(pvalue = 1)@result, "richResult")
  expect_gt(nrow(l@result), 0)

  expect_warning(tt <- testedTerms(l), "tested table unavailable: .*no tested table was recorded")
  expect_null(tt)
  expect_identical(enrichmentInfo(l), unavailable_record)
  expect_identical(result(l), as.data.frame(l@result))
  expect_identical(compareResult(list(a = l)), compareResult(list(a = l@result)))
})

test_that("a result stripped of its tested slot and metadata is unavailable, not reconstructed", {
  x <- .fora_reference(pvalue = 1)
  attr(x, "tested") <- NULL
  attr(x, "richR_metadata") <- NULL
  expect_false(methods::.hasSlot(x, "tested"))
  expect_gt(nrow(x@result), 0)

  expect_warning(tt <- testedTerms(x), "tested table unavailable: .*no tested slot")
  expect_null(tt)
  expect_identical(enrichmentInfo(x), unavailable_record)
  expect_identical(result(x), as.data.frame(x@result))
  expect_identical(compareResult(list(a = x)), compareResult(list(a = x@result)))
})
