# Shared fixtures and base-R readers of the report tests (SPEC-RICHR-008),
# built from base R, testthat and the public API, with no network and no
# package outside DESCRIPTION (report.html is read with regular expressions,
# which is exact because richReport() double-quotes every attribute value and
# escapes every text '<'; test-report.R checks that grammar).
#
# A80    annotation over G01..G80: K (G01-G08, "kinase signalling"), M
#        (G01-G06, G09), L (G11-G16), N (G30-G37), W (G17-G29), Z (G41-G80);
#        labels mix upper and lower case, so the session collation changes
#        the order of a figure's term axis (plan ledger L-S12).
# F-Q    q1 reports K, M, L (W tested only, N untested); q2 reports N (K
#        tested only, M, L, W untested).
# F-QM   q1's query mapped by mapGenes(), with an unmapped and a duplicate id.
# F-QX   a mapped query whose mapping record no longer matches the input.
# F-LEG  as(q1@result, "richResult"): no analysis record.
# F-NEW  new("richResult", ...): refused by writeRich().
# F-FIG  12 reported sets in result order S01..S12, 4 with Pvalue > 0.05.
# F-ZERO a run that reports nothing.
# F-SENS a 2 x 2 richSensitivity() grid over F-Q's query.
#
# Every builder is paired with a role check; tests assert the roles first, so
# a fixture drift fails loudly instead of silently changing what is checked.

.rp_G <- function(i) sprintf("G%02d", i)

.rp_a80 <- function(rename = NULL) {
  G <- .rp_G(1:80)
  if (!is.null(rename)) G[G == names(rename)] <- rename[[1]]
  data.frame(GeneID = c(G[1:8], G[c(1:6, 9)], G[11:16], G[30:37], G[17:29], G[41:80]),
             Term = rep(c("K", "M", "L", "N", "W", "Z"), c(8, 7, 6, 8, 13, 40)),
             Annot = rep(c("kinase signalling", "M set", "L set", "N set", "W set",
                           "Z set"), c(8, 7, 6, 8, 13, 40)),
             stringsAsFactors = FALSE)
}

.rp_enrich <- function(q, annot = .rp_a80(), pvalue = 0.05) {
  suppressMessages(richR::enrich(q, annot, pvalue = pvalue, padj = NULL, minSize = 2,
                                 minGSSize = 3))
}

.rp_q1_query <- function() .rp_G(c(1:6, 11:14, 17))

.rp_fq <- function() {
  list(q1 = .rp_enrich(.rp_q1_query()), q2 = .rp_enrich(.rp_G(c(7, 8, 30:35))))
}

.rp_fq_roles <- function(q1, q2) {
  b <- richR::bindRich(list(q1 = q1, q2 = q2))
  st <- function(id, g) b$state[b$Annot == id & as.character(b$group) == g]
  testthat::expect_identical(
    c(st("K", "q1"), st("M", "q1"), st("L", "q1"), st("W", "q1"), st("N", "q1")),
    c("reported", "reported", "reported", "tested_not_reported", "not_tested"))
  testthat::expect_identical(
    c(st("N", "q2"), st("K", "q2"), st("M", "q2"), st("L", "q2"), st("W", "q2")),
    c("reported", "tested_not_reported", "not_tested", "not_tested", "not_tested"))
}

.rp_map_table <- function() {
  G <- .rp_G(1:80)
  data.frame(from = tolower(G), to = G, stringsAsFactors = FALSE)
}

.rp_fqm <- function() {
  v <- suppressMessages(richR::mapGenes(c(tolower(.rp_q1_query()), "zzz", "g01"),
                                        .rp_map_table(), "from", "to"))
  .rp_enrich(v)
}

.rp_fqx <- function() {
  v <- suppressMessages(richR::mapGenes(tolower(.rp_q1_query()), .rp_map_table(),
                                        "from", "to"))
  v[11] <- "G18"
  testthat::expect_warning(qx <- .rp_enrich(v), "does not match the submitted input")
  qx
}

.rp_fleg <- function(q1) methods::as(q1@result, "richResult")

.rp_fnew <- function(q1) methods::new("richResult", result = q1@result, sep = ",")

.rp_ffig <- function() {
  H <- sprintf("H%03d", 1:1000)
  fa <- do.call(rbind, lapply(1:12, function(i) {
    data.frame(GeneID = H[((i - 1) * 20 + 1):(i * 20)], Term = sprintf("S%02d", i),
               Annot = sprintf("set %02d", i), stringsAsFactors = FALSE)
  }))
  fa <- rbind(fa, data.frame(GeneID = H[241:1000], Term = "BG", Annot = "background set",
                             stringsAsFactors = FALSE))
  q <- unlist(lapply(1:12, function(i) H[((i - 1) * 20 + 1):((i - 1) * 20 + 13 - i)]))
  suppressMessages(richR::enrich(q, fa, pvalue = 1, padj = NULL, minSize = 1,
                                 minGSSize = 3, maxGSSize = 100))
}

.rp_fzero <- function() .rp_enrich(.rp_G(41:42), pvalue = 1e-10)

.rp_fsens <- function() {
  suppressMessages(richR::richSensitivity(
    .rp_q1_query(), .rp_a80(), backgrounds = list(all = NULL, part = .rp_G(1:60)),
    maxGSSize = c(10, 50), minSize = 2, minGSSize = 3, pvalue = 0.05, padj = NULL))
}

# TH: themes of R-MAIN with the first label edited; ST: statistics for q1.
.rp_th <- function(q1, q2) {
  th <- richR::richThemes(list(q1 = q1, q2 = q2, f = NULL))
  th$themes$label[1] <- "Edited theme"
  th
}

.rp_st <- function() {
  list(q1 = data.frame(gene = .rp_G(1:7), log2FC = seq(-1.5, 1.5, length.out = 7),
                       stringsAsFactors = FALSE))
}

# R-MAIN: the report every AC-001 / AC-002 / AC-008 leg reads.
.rp_main <- function(path, q = .rp_fq(), ss = .rp_fsens()) {
  th <- .rp_th(q$q1, q$q2)
  richR::richReport(list(q1 = q$q1, q2 = q$q2, f = NULL), path, themes = th,
                    stats = .rp_st(), sensitivity = ss)
  list(q = q, th = th, ss = ss)
}

# A fresh scratch folder, removed when the calling test ends.
.rp_tmp <- function(env = parent.frame()) {
  d <- tempfile("report-test-")
  dir.create(d)
  do.call(on.exit, list(substitute(unlink(d, recursive = TRUE), list(d = d)), add = TRUE),
          envir = env)
  normalizePath(d)
}

# ---- base-R readers -------------------------------------------------------

.rp_html <- function(p) {
  paste(readLines(file.path(p, "report.html"), warn = FALSE, encoding = "UTF-8"),
        collapse = "\n")
}

.rp_unescape <- function(s) {
  s <- gsub("&lt;", "<", s, fixed = TRUE)
  s <- gsub("&gt;", ">", s, fixed = TRUE)
  s <- gsub("&quot;", "\"", s, fixed = TRUE)
  s <- gsub("&#39;", "'", s, fixed = TRUE)
  gsub("&amp;", "&", s, fixed = TRUE)
}

# Visible text: tags dropped, entities unescaped, whitespace runs collapsed.
.rp_text <- function(html) {
  s <- .rp_unescape(gsub("<[^>]*>", " ", html))
  gsub("\\s+", " ", s)
}

.rp_has <- function(html, s) grepl(s, .rp_text(html), fixed = TRUE)

.rp_match_all <- function(x, re) regmatches(x, gregexpr(re, x, perl = TRUE))[[1]]

# The HTML from the element with id `id` up to the next h2 (or the end).
.rp_section <- function(html, id, until = "<h2") {
  start <- regexpr(sprintf(" id=\"%s\"", id), html, fixed = TRUE)
  if (start < 0) return(NA_character_)
  rest <- substring(html, start + 1L)
  stop_at <- regexpr(until, rest, fixed = TRUE)
  if (stop_at > 0) substring(rest, 1L, stop_at - 1L) else rest
}

# Rows of the first table after `marker` (text) in `html`, as lists of cell
# texts; row 1 is the header.
.rp_table_rows <- function(html, marker = NULL) {
  if (!is.null(marker)) {
    pos <- regexpr(marker, html, fixed = TRUE)
    if (pos < 0) return(NULL)
    html <- substring(html, pos)
  }
  tab <- regmatches(html, regexpr("(?s)<table>.*?</table>", html, perl = TRUE))
  if (!length(tab)) return(NULL)
  rows <- .rp_match_all(tab, "(?s)<tr[^>]*>.*?</tr>")
  lapply(rows, function(r) {
    cells <- .rp_match_all(r, "(?s)<t[dh]>.*?</t[dh]>")
    trimws(.rp_text(cells))
  })
}

# The value cell of the record row whose key cell is `key`.
.rp_record_value <- function(html, key) {
  rows <- .rp_table_rows(html)
  hit <- Filter(function(r) length(r) == 2L && identical(r[1], key), rows)
  if (length(hit) != 1L) return(NA_character_)
  hit[[1]][2]
}

.rp_captions <- function(html) {
  trimws(.rp_text(.rp_match_all(html, "<figcaption>.*?</figcaption>")))
}

.rp_ids <- function(html) {
  sub("^ id=\"(.*)\"$", "\\1", .rp_match_all(html, " id=\"[^\"]*\""))
}

.rp_hrefs <- function(html) {
  sub("^href=\"(.*)\"$", "\\1", .rp_match_all(html, "href=\"[^\"]*\""))
}

# The report's own anchor rule, computed independently: g<i>- and the
# lowercase hex of the identifier's UTF-8 bytes.
.rp_anchor <- function(i, id) {
  paste0("g", i, "-", vapply(id, function(v) {
    paste(as.character(charToRaw(enc2utf8(v))), collapse = "")
  }, character(1), USE.NAMES = FALSE))
}

.rp_manifest <- function(p) {
  jsonlite::fromJSON(paste(readLines(file.path(p, "manifest.json"), warn = FALSE,
                                     encoding = "UTF-8"), collapse = "\n"),
                     simplifyVector = FALSE)
}

# Decode a report table with the record manifest.json holds for it, through
# the reader richReport() itself uses.
.rp_decode <- function(p, stem) richR:::.report_read_table(p, stem)

.rp_files <- function(p) {
  sort(list.files(p, recursive = TRUE, all.files = TRUE), method = "radix")
}

.rp_md5 <- function(p) tools::md5sum(file.path(p, .rp_files(p)))

.rp_copy <- function(from, to) {
  dir.create(to)
  file.copy(list.files(from, full.names = TRUE, all.files = TRUE, no.. = TRUE), to,
            recursive = TRUE)
  invisible(to)
}

.rp_rscript <- function() file.path(R.home("bin"), "Rscript")

# AC-008: q holds p's file names; tables/, results/ and rerun.R are byte
# identical; report.html and manifest.json are identical once p's generation
# time (which occurs exactly once in each) is replaced by q's.
.rp_expect_rerun_identical <- function(p, q) {
  files <- .rp_files(p)
  testthat::expect_identical(.rp_files(q), files)
  fixed <- grep("^(tables/|results/|rerun\\.R$)", files, value = TRUE)
  testthat::expect_gt(length(fixed), 10)
  testthat::expect_identical(unname(tools::md5sum(file.path(q, fixed))),
                             unname(tools::md5sum(file.path(p, fixed))))
  gp <- .rp_manifest(p)$generated
  gq <- .rp_manifest(q)$generated
  stamp <- "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$"
  testthat::expect_true(grepl(stamp, gp) && grepl(stamp, gq))
  for (f in c("report.html", "manifest.json")) {
    lp <- readLines(file.path(p, f), warn = FALSE, encoding = "UTF-8")
    hits <- regmatches(lp, gregexpr(gp, lp, fixed = TRUE))
    testthat::expect_identical(sum(lengths(hits)), 1L, label = f)
    testthat::expect_identical(sub(gp, gq, lp, fixed = TRUE),
                               readLines(file.path(q, f), warn = FALSE, encoding = "UTF-8"),
                               label = f)
  }
}

.rp_lib_env <- function() {
  paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
}
