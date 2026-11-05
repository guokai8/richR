# Live KEGG pathway labels (builtin = FALSE) must survive an unreachable or
# incomplete KEGG listing: fall back to the bundled labels and say so.
# Catches: a crash when KEGG is unreachable, a silent fallback, and a spurious
# warning on a healthy listing. keggList is mocked; nothing reaches the network.

collect_warnings <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = messages)
}

kegg_listing_cases <- list(
  failing = function(...) stop("simulated KEGG outage"),
  without_04148 = function(...) c(map00010 = "live Glycolysis"),
  healthy = function(...) c(map00010 = "live Glycolysis", map04148 = "live Efferocytosis")
)

test_that("live KEGG labels fall back to the bundled labels with a warning naming the reason", {
  bundled <- richR:::.get_kg_dat(builtin = TRUE)
  expect_true("04148" %in% rownames(bundled))

  failing <- collect_warnings(testthat::with_mocked_bindings(
    keggList = kegg_listing_cases$failing,
    richR:::.get_kg_dat(builtin = FALSE),
    .package = "richR"
  ))
  expect_identical(failing$value, bundled)
  expect_equal(length(failing$warnings), 1)
  expect_true(grepl("simulated KEGG outage", failing$warnings[1], fixed = TRUE))
  expect_true(grepl("bundled", failing$warnings[1], fixed = TRUE))

  incomplete <- collect_warnings(testthat::with_mocked_bindings(
    keggList = kegg_listing_cases$without_04148,
    richR:::.get_kg_dat(builtin = FALSE),
    .package = "richR"
  ))
  expect_identical(incomplete$value, bundled)
  expect_equal(length(incomplete$warnings), 1)
  expect_true(grepl("04148", incomplete$warnings[1], fixed = TRUE))
  expect_true(grepl("bundled", incomplete$warnings[1], fixed = TRUE))

  healthy <- collect_warnings(testthat::with_mocked_bindings(
    keggList = kegg_listing_cases$healthy,
    richR:::.get_kg_dat(builtin = FALSE),
    .package = "richR"
  ))
  expect_equal(length(healthy$warnings), 0)
  expect_equal(healthy$value[c("00010", "04148"), "annotation"],
               c("live Glycolysis", "live Efferocytosis"))
})

test_that("richKEGG with builtin = FALSE labels from the bundled snapshot when KEGG is unreachable", {
  bundled <- richR:::.get_kg_dat(builtin = TRUE)
  annotation <- data.frame(
    GeneID = c("A", "B", "C", "D"),
    PATH = c("00010", "00010", "04148", "04148"),
    stringsAsFactors = FALSE
  )
  run <- function(listing) {
    collect_warnings(testthat::with_mocked_bindings(
      keggList = listing,
      richR::richKEGG(c("A", "C"), annotation, pvalue = 1, padj = 1,
                      minSize = 1, minGSSize = 1, builtin = FALSE),
      .package = "richR"
    ))
  }
  label_of <- function(result, id) {
    as.character(result@result$Term[result@result$Annot == id])
  }

  for (case in c("failing", "without_04148")) {
    fallback <- run(kegg_listing_cases[[case]])
    expect_s4_class(fallback$value, "richResult")
    expect_equal(length(fallback$warnings), 1)
    expect_true(grepl("bundled", fallback$warnings[1], fixed = TRUE))
    expect_equal(label_of(fallback$value, "00010"), bundled["00010", "annotation"])
    expect_equal(label_of(fallback$value, "04148"), bundled["04148", "annotation"])
  }

  healthy <- run(kegg_listing_cases$healthy)
  expect_equal(length(healthy$warnings), 0)
  expect_equal(label_of(healthy$value, "00010"), "live Glycolysis")
  expect_equal(label_of(healthy$value, "04148"), "live Efferocytosis")
})
