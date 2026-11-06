# Shared fixtures of the tested-table tests (SPEC-RICHR-004), built with the
# public API and no network.
#
# F-ORA: a 20-gene universe and five gene sets in four roles under the
# reference call (.fora_reference(): pvalue = 0.05, padj = NULL,
# minSize = 2, minGSSize = 3, keepRich = FALSE):
#   T1 reported (overlap 5 of 5);
#   T2 tested, overlap 2 >= minSize, fails the p-value cutoff only;
#   T5 tested, overlap 1 < minSize;
#   T3 size-eligible, overlap 0;
#   T4 size 2 < minGSSize, never tested.
# fixtures/make-legacy-fixtures.R builds the same annotation for the legacy
# fixtures. Tests assert each role before relying on it (.fora_roles()).

.fora_annot <- function() {
  genes <- sprintf("G%02d", 1:20)
  mk <- function(term, label, idx) {
    data.frame(GeneID = genes[idx], Term = term, Annot = label,
               stringsAsFactors = FALSE)
  }
  rbind(
    mk("T1", "reported set", 1:5),
    mk("T2", "tested set failing the cutoff", c(6:10, 12:14)),
    mk("T3", "eligible set without overlap", 18:20),
    mk("T4", "set below minGSSize", c(1, 15)),
    mk("T5", "tested set below minSize", c(11, 15:17))
  )
}

.fora_query <- function() sprintf("G%02d", c(1:7, 11))

.fora_reference <- function(...) {
  args <- utils::modifyList(
    list(pvalue = 0.05, padj = NULL, minSize = 2, minGSSize = 3,
         keepRich = FALSE),
    list(...)
  )
  suppressMessages(do.call(
    richR::enrich, c(list(.fora_query(), .fora_annot()), args)
  ))
}

# Asserts the F-ORA roles on a tested table so a fixture drift fails loudly.
.fora_roles <- function(tested, pvalue = 0.05) {
  testthat::expect_setequal(tested$Annot, c("T1", "T2", "T4", "T5"))
  sig <- stats::setNames(tested$Significant, tested$Annot)
  p <- stats::setNames(tested$Pvalue, tested$Annot)
  testthat::expect_true(sig[["T1"]] >= 2 && p[["T1"]] <= pvalue)
  testthat::expect_true(sig[["T2"]] >= 2 && p[["T2"]] > pvalue)
  testthat::expect_identical(as.numeric(sig[["T5"]]), 1)
}

# F-GSEA: the characterization GSEA fixtures under a fixed seed.
.fgsea_ranked <- function(name = "gsea_ranked_a.csv") {
  df <- .char_read(name)
  stats::setNames(df$score, df$GeneID)
}

.fgsea_run <- function(name = "gsea_ranked_a.csv", ...,
                       annot = .char_read("gsea_annot.csv", character = TRUE)) {
  out <- NULL
  utils::capture.output(
    out <- richR::richGSEA(.fgsea_ranked(name), annot[, 1:2, drop = FALSE], seed = 42, nproc = 1,
                           ...)
  )
  if (identical(list(...)$pvalue, 1)) {
    # This fixture intentionally displays all backend-returned rows. Public
    # upstream strict cutoffs are tested separately; no statistics change.
    full <- out@tested
    out@result <- full[, setdiff(names(full), c("reported", "Term")), drop = FALSE]
    out@tested$reported <- TRUE
  }
  out
}
