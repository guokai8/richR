stage_mapping_vector <- function(decisions) {
  kept <- decisions$decision == "mapped"
  value <- if (is.null(decisions$score)) decisions$target[kept] else
    stats::setNames(decisions$score[kept], decisions$target[kept])
  attr(value, "richR_mapping") <- list(decisions = decisions,
                                      source = list(kind = "table"))
  value
}

# Captured from installed upstream dca5857, schema 1, before occurrence capture.
stage_legacy_ora <- function() {
  readRDS(testthat::test_path("fixtures", "legacy-audit-plain-v012.rds"))
}

test_that("audit frames preserve column types, list candidates and zero rows", {
  decisions <- list(position = 1:3, input_id = c("a", "b", "c"),
                    candidates = list(c("A", "A2"), character(), "C"),
                    target = c("A", NA, "C"),
                    decision = c("mapped", "unmapped", "mapped"),
                    reason = c("first candidate", "no candidate", "one candidate"))
  ranking <- decisions
  ranking$score <- c(2, NA_real_, -1)
  empty <- lapply(decisions, function(column) column[FALSE])
  for (record in list(decisions, ranking, empty)) {
    value <- stage_mapping_vector(record)
    before <- serialize(value, NULL)
    audit <- mappingAudit(value)
    expected <- data.frame(position = record$position, input_id = record$input_id,
                           target = record$target, decision = record$decision,
                           reason = record$reason, stringsAsFactors = FALSE)
    expected$candidates <- record$candidates
    if (!is.null(record$score)) expected$score <- record$score
    expected <- expected[names(record)]
    expect_identical(audit, expected)
    # R 3.5 has no list2DF; compare with it only on runtimes that provide it.
    canonical <- get0("list2DF", envir = baseenv(), inherits = FALSE)
    if (is.function(canonical)) expect_identical(audit, canonical(record))
    expect_identical(names(audit), names(record))
    expect_true(is.list(audit$candidates))
    expect_identical(.row_names_info(audit), -as.integer(length(record$target)))
    expect_identical(serialize(value, NULL), before)
  }
})

test_that("conversion does not change either the input or mapping source", {
  mapping <- data.frame(source = c("a", "b", "b", "c"),
                        target = c("A", "B1", "B2", "A"))
  ids <- c("a", "b", "c", NA_character_)
  before <- serialize(list(ids, mapping), NULL)
  value <- mapGenes(ids, mapping, "source", "target")
  expect_identical(as.vector(value), "A")
  expect_identical(mappingAudit(value)$decision,
                   c("mapped", "ambiguous", "merged", "invalid"))
  expect_identical(serialize(list(ids, mapping), NULL), before)
})

test_that("unavailable and malformed vector records explain the missing evidence", {
  expect_warning(value <- mappingAudit(c("a", "b")), "no mapping record")
  expect_null(value)
  malformed <- list(
    1L, list(), list(decisions = list(target = "A", decision = 1L)),
    list(decisions = list(target = c("A", "B"), decision = "mapped")),
    list(decisions = list(decision = "mapped")))
  for (record in malformed) {
    object <- structure("A", richR_mapping = record)
    expect_warning(value <- mappingAudit(object), "mapping record is malformed")
    expect_null(value)
  }
})

test_that("legacy ORA audit reads recorded facts without inventing occurrences", {
  object <- stage_legacy_ora()
  before <- serialize(object, NULL)
  expect_no_warning(audit <- mappingAudit(object))
  query <- audit[audit$role == "query", , drop = FALSE]
  background <- audit[audit$role == "background", , drop = FALSE]
  expect_identical(query$id, c("A", "C", "X", NA, NA))
  expect_identical(query$decision,
                   c("effective", "effective", "not_annotated", "invalid", "invalid"))
  expect_true(all(is.na(query$occurrences[!is.na(query$id)])))
  expect_identical(query$occurrences[is.na(query$id)], c(1L, 1L))
  expect_true(all(is.na(query$mapping)))
  expect_identical(background$id, c("A", "B", "C", "D", "X"))
  expect_identical(background$decision,
                   c(rep("effective", 4), "not_annotated"))
  expect_true(all(is.na(background$occurrences)))
  expect_identical(serialize(object, NULL), before)

  # The audit is evidence from metadata, independently of displayed tables.
  altered <- object
  altered@result <- altered@result[FALSE, , drop = FALSE]
  altered@gene <- "fabricated"
  altered@detail <- data.frame(GeneID = "fabricated")
  expect_identical(mappingAudit(altered), audit)
})

test_that("legacy GSEA audit retains recorded ranking, scores and ties", {
  ranking <- stats::setNames(c(3, 2, 2, 1, .5, 0, -.5, -1, -2, -3, -4, -5),
                             paste0("G", 1:12))
  annotation <- data.frame(GeneID = paste0("G", 1:10),
                           Term = rep(c("T1", "T2"), each = 5))
  # Genuine pre-capture record from installed upstream dca5857, schema 1.
  object <- readRDS(testthat::test_path("fixtures", "legacy-audit-gsea-v012.rds"))
  before <- serialize(object, NULL)
  expect_no_warning(audit <- mappingAudit(object))
  expect_identical(audit$id, names(ranking))
  expect_identical(audit$score, unname(ranking))
  expect_identical(audit$tied, c(FALSE, TRUE, TRUE, rep(FALSE, 9)))
  expect_identical(audit$decision,
                   c(rep("ranked_annotated", 10), rep("ranked_not_annotated", 2)))
  expect_true(all(audit$role == "ranked"))
  expect_true(all(is.na(audit$mapping)))
  expect_true(all(is.na(audit$occurrences)))
  expect_identical(serialize(object, NULL), before)
})

test_that("missing or unavailable result metadata gives a reason rather than guesses", {
  object <- stage_legacy_ora()
  attr(object, "richR_metadata") <- NULL
  expect_warning(value <- mappingAudit(object), "metadata not recorded")
  expect_null(value)
  attr(object, "richR_metadata") <- list(available = FALSE,
                                         reason = "record unavailable for fixture")
  expect_warning(value <- mappingAudit(object), "record unavailable for fixture")
  expect_null(value)
})
