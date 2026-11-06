# writeRich() / readRich() (SPEC-RICHR-004 REQ-011..REQ-014, AC-005 I/O leg,
# AC-010..AC-014): a result object stored as a folder of TSV tables plus a JSON
# manifest. Fixtures are built with the public API (helper-tested.R) and every
# folder is written under tempdir().

io_dir <- function() tempfile("rich-io-")

# F-S22 (AC-010): 22 gene sets with five-digit identifiers holding leading
# zeros, 4 genes each in a 100-gene universe; the query overlaps 18 sets in 2
# genes and 4 sets in 1 gene, so enrich() reports 18 rows and records 4
# tested-but-not-reported rows.
fs22 <- function() {
  u <- sprintf("g%03d", 1:100)
  annot <- do.call(rbind, lapply(1:22, function(i) {
    data.frame(GeneID = u[(4 * i - 3):(4 * i)], Term = sprintf("%05d", 10 * i),
               Annot = paste("set", i), stringsAsFactors = FALSE)
  }))
  query <- c(unlist(lapply(1:18, function(i) u[(4 * i - 3):(4 * i - 2)])),
             u[4 * (19:22) - 3])
  suppressMessages(richR::enrich(query, annot, universe = u, pvalue = 1,
                                 padj = NULL, minSize = 2, minGSSize = 1))
}

# The 17 character values of acceptance.md L-P4.
lp4_values <- c(NA, "", "NA", "a\tb", "l1\nl2", "c\rr", "q\"uote", "it's",
                "back\\slash", "\\N", "\\\\N", "#hash", " lead", "trail ",
                "日本", "00010", "\\t")

# F-S22 edited by slot assignment only: special doubles, the L-P4 labels in
# @result and the same labels in @tested for the same identifiers, and a
# factor column in @detail.
fs22_edited <- function() {
  x <- fs22()
  x@result$zscore[1:6] <- c(NA, NaN, Inf, -Inf, 5e-324, 1 / 3)
  x@result$Term[1:17] <- lp4_values
  x@tested$Term[match(x@result$Annot[1:17], x@tested$Annot)] <- lp4_values
  x@detail$flag <- factor(rep_len(c("b", "a", NA), nrow(x@detail)),
                          levels = c("b", "a", "z"))
  x
}

kegg_result <- function() {
  suppressMessages(richKEGG(
    .char_read("ora_query.csv")$GeneID,
    .char_read("kegg_annot.csv", character = TRUE),
    pvalue = 0.05, padj = NULL, minSize = 2, minGSSize = 3
  ))
}

read_manifest <- function(p) {
  jsonlite::fromJSON(file.path(p, "manifest.json"), simplifyVector = FALSE)
}

dir_state <- function(p) {
  files <- sort(list.files(p, all.files = TRUE, no.. = TRUE))
  stats::setNames(unname(tools::md5sum(file.path(p, files))), files)
}

# Round-trip equality of acceptance.md § D.0: per table identical column
# names, classes, row names and columns (doubles bit-exact); a zero-column
# table identical as a whole; every other slot identical; the tested table
# unavailable exactly when it was; enrichmentInfo() identical.
expect_same_table <- function(b, a, what) {
  expect_identical(names(b), names(a), info = what)
  expect_identical(lapply(b, class), lapply(a, class), info = what)
  expect_identical(rownames(b), rownames(a), info = what)
  for (j in seq_along(a)) {
    expect_identical(b[[j]], a[[j]], info = paste(what, names(a)[j]))
  }
  if (!ncol(a)) expect_identical(b, a, info = what)
}

tested_available <- function(x) !is.null(suppressWarnings(testedTerms(x)))

expect_round_trip <- function(b, a) {
  expect_identical(class(b), class(a))
  tables <- if (methods::is(a, "richResult")) c("result", "detail") else "result"
  for (s in tables) expect_same_table(methods::slot(b, s), methods::slot(a, s), s)
  expect_identical(tested_available(b), tested_available(a))
  if (tested_available(a)) expect_same_table(b@tested, a@tested, "tested")
  for (s in setdiff(methods::slotNames(class(a)), c(tables, "tested"))) {
    expect_identical(methods::slot(b, s), methods::slot(a, s), info = s)
  }
  expect_identical(enrichmentInfo(b), enrichmentInfo(a))
}

round_trip <- function(x) {
  p <- io_dir()
  writeRich(x, p)
  readRich(p)
}

copy_folder <- function(src) {
  d <- io_dir()
  dir.create(d)
  file.copy(list.files(src, full.names = TRUE), d)
  d
}

write_manifest <- function(p, m) {
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, null = "null",
                              digits = I(17), pretty = TRUE),
             file.path(p, "manifest.json"))
}

edit_manifest <- function(p, f) write_manifest(p, f(read_manifest(p)))

# Edits the typed metadata the manifest embeds (jsonlite::serializeJSON form).
edit_metadata <- function(p, f) {
  edit_manifest(p, function(m) {
    meta <- jsonlite::unserializeJSON(jsonlite::toJSON(
      m$metadata, auto_unbox = TRUE, null = "null", digits = I(17)))
    m$metadata <- jsonlite::fromJSON(
      jsonlite::serializeJSON(f(meta), digits = I(17)), simplifyVector = FALSE)
    m
  })
}

edit_lines <- function(p, file, f) {
  path <- file.path(p, file)
  writeLines(f(readLines(path, encoding = "UTF-8")), path, useBytes = TRUE)
  # Keep checksum valid here so each test reaches its intended typed parser.
  edit_manifest(p, function(m) {
    m$tables <- lapply(m$tables, function(e) {
      if (identical(e$file, file)) e$sha256 <- digest::digest(path, algo="sha256", file=TRUE)
      e
    })
    m
  })
}

# ---- writeRich: folder content (AC-010, AC-012) -----------------------------

test_that("writeRich writes the manifest and identifiers verbatim", {
  # Catches identifiers converted to numbers (L-P2: "00010" read as 10) and a
  # manifest that does not say what the folder holds.
  x <- fs22()
  expect_identical(nrow(x@result), 18L)
  expect_identical(sum(!testedTerms(x)$reported), 4L)
  expect_true(all(grepl("^[0-9]{5}$", x@result$Annot)))
  p <- io_dir()
  expect_invisible(writeRich(x, p))

  expect_setequal(list.files(p), c("manifest.json", "result.tsv",
                                   "tested.tsv", "detail.tsv"))
  m <- read_manifest(p)
  expect_identical(m$format, "richR-result")
  expect_identical(m$format_version, 2L)
  expect_identical(m$class, "richResult")
  expect_true(m$tested_recorded)
  ids <- utils::read.delim(file.path(p, "result.tsv"), quote = "",
                           na.strings = character(0), comment.char = "",
                           colClasses = "character")$Annot
  expect_identical(ids, x@result$Annot)
})

test_that("writeRich records a zero-column table in the manifest only", {
  # Catches a header-only file for a zero-column table, which cannot be read
  # back (L-P5).
  x <- .fora_reference()
  x@detail <- x@detail[, 0, drop = FALSE]
  expect_gt(nrow(x@detail), 0)
  p <- io_dir()
  writeRich(x, p)
  expect_false(file.exists(file.path(p, "detail.tsv")))
  m <- read_manifest(p)
  detail <- Filter(function(t) identical(t$slot, "detail"), m$tables)[[1]]
  expect_null(detail$file)
  expect_identical(detail$ncol, 0L)
  expect_identical(detail$nrow, nrow(x@detail))
})

# ---- writeRich: refusals before any file operation (AC-005, AC-014) --------

test_that("writeRich refuses a grouped or incomplete view and creates nothing", {
  # Catches a grouped view silently stored ungrouped, and a view without its
  # identifier or label column stored as if complete (REQ-012).
  r <- .fora_reference()
  p <- io_dir()
  expect_error(writeRich({ z <- r; z@result <- dplyr::group_by(r@result, Term); z }, p), "@result.*grouped_df")
  expect_false(file.exists(p))
  expect_error(writeRich(dplyr::select(r, -Term), p), "Term")
  expect_false(file.exists(p))
  expect_error(writeRich(dplyr::select(r, -Annot), p), "Annot")
  expect_false(file.exists(p))
  g <- .fgsea_run(minSize = 3, maxSize = 9, pvalue = 0.05)
  expect_error(writeRich(dplyr::select(g, -pathway), p), "pathway")
  expect_false(file.exists(p))
})

test_that("writeRich never touches a folder it did not write", {
  # Catches user files deleted or mixed into a result folder, a silent
  # overwrite of an earlier result, and an unsupported column silently
  # coerced to text (REQ-012).
  x <- .fora_reference()
  d <- io_dir()
  dir.create(d)
  writeLines("user notes", file.path(d, "keep.txt"))
  before <- dir_state(d)
  expect_error(writeRich(x, d), "not empty")
  expect_identical(dir_state(d), before)
  expect_error(writeRich(x, d, overwrite = TRUE), "not written by writeRich")
  expect_identical(dir_state(d), before)

  w <- io_dir()
  writeRich(x, w)
  written <- dir_state(w)
  expect_error(writeRich(.fora_reference(pvalue = 1), w), "not empty")
  expect_identical(dir_state(w), written)

  y <- x
  y@result$when <- as.Date("2026-11-05")
  e <- io_dir()
  expect_error(writeRich(y, e), "when")
  expect_false(file.exists(e))
  y <- x
  y@tested$score <- 1:nrow(y@tested) + 0i
  expect_error(writeRich(y, e), "score")
  expect_false(file.exists(e))
  expect_error(writeRich(x@result, e), "richResult or GSEAResult")
  expect_false(file.exists(e))
})

# ---- readRich: round trip (AC-005, AC-010..AC-012, AC-014) ------------------

test_that("readRich restores an ORA result exactly, statistics and labels included", {
  # Catches precision loss in stored statistics (L-P1), identifiers converted
  # to numbers (L-P2), NA / "" / "NA" conflated or rows lost (L-P4b), a TAB or
  # newline splitting a cell, factor levels lost, extra columns of a pipeline
  # result dropped, and metadata dropped or retyped on reload (L-P3).
  x <- fs22_edited()
  expect_no_error(methods::validObject(x))
  b <- round_trip(x)
  expect_round_trip(b, x)
  expect_identical(b@result$zscore[1:6], c(NA, NaN, Inf, -Inf, 5e-324, 1 / 3))
  expect_identical(b@result$Term[1:17], lp4_values)
  expect_identical(levels(b@detail$flag), c("b", "a", "z"))

  k <- kegg_result()
  expect_true(all(c("ko", "Level3", "Level2", "Level1") %in% names(k@result)))
  expect_round_trip(round_trip(k), k)
})

test_that("readRich restores GSEA results with their ranking, leading edges and provenance", {
  # Catches a GSEA ranking reordered or renamed on reload, a list column
  # flattened or split wrongly, an empty leading edge reloaded as list(), and
  # GSEA provenance dropped.
  skip_if_not_installed("fgsea")
  ranked <- .fgsea_ranked()
  ranked[["G05"]] <- ranked[["G04"]]
  run <- function(annot, ...) {
    out <- NULL
    utils::capture.output(out <- suppressWarnings(richR::richGSEA(
      ranked, annot, seed = 42, nproc = 1, minSize = 3, maxSize = 9,
      pvalue = 0.05, ...)))
    out
  }
  annot <- .char_read("gsea_annot.csv", character = TRUE)
  gt <- run(annot, table = TRUE)
  expect_gt(anyDuplicated(unname(gt@input)), 0)
  expect_true(any(gt@input < 0))
  gl <- run(annot, table = FALSE)
  expect_gt(nrow(gl@result), 0)
  expect_type(gl@result$leadingEdge, "list")
  gl@result$leadingEdge[[1]] <- character(0)
  ak <- annot
  ak$Term <- unname(c(T1 = "00010", T2 = "00061", T3 = "00190", T4 = "00230",
                      T5 = "00970")[ak$Term])
  gk <- run(ak, KEGG = TRUE)
  expect_true("Level1" %in% names(testedTerms(gk)))

  for (g in list(gt, gl, gk)) expect_round_trip(round_trip(g), g)
  bl <- round_trip(gl)
  expect_identical(bl@result$leadingEdge[[1]], character(0))
  expect_identical(names(bl@input), names(gl@input))
  info <- enrichmentInfo(bl)$provenance
  expect_false(is.null(info$fgsea_version))
  expect_false(is.null(info$rng))
})

test_that("readRich restores empty results, zero-column tables and unavailable records", {
  # Catches an empty result rejected or reloaded with logical columns, a
  # zero-column table that cannot be read back (L-P5) or loses its row count,
  # and a legacy object reloaded with a fabricated tested table or upgraded
  # metadata.
  none <- .fora_reference(pvalue = 1e-12)
  expect_identical(nrow(none@result), 0L)
  expect_gt(nrow(testedTerms(none)), 0)
  expect_round_trip(round_trip(none), none)

  emptied <- .fora_reference()
  emptied@tested <- emptied@tested[0, ]
  expect_round_trip(round_trip(emptied), emptied)

  zero <- .fora_reference()
  zero@detail <- zero@detail[, 0, drop = FALSE]
  expect_round_trip(round_trip(zero), zero)
  named <- .fora_reference()
  named@detail <- data.frame(a = 1:3, row.names = c("x", "y", "z"))[, 0, drop = FALSE]
  expect_round_trip(round_trip(named), named)

  as_obj <- methods::as(.fora_reference()@result, "richResult")
  expect_identical(dim(as_obj@detail), c(0L, 0L))
  p <- io_dir()
  writeRich(as_obj, p)
  expect_false(file.exists(file.path(p, "detail.tsv")))
  b <- readRich(p)
  expect_round_trip(b, as_obj)
  expect_false(enrichmentInfo(b)$available)

  for (name in c("legacy_richResult_fca5276.rds", "legacy_GSEAResult_fca5276.rds")) {
    legacy <- readRDS(test_path("fixtures", name))
    p <- io_dir()
    writeRich(legacy, p)
    expect_false(file.exists(file.path(p, "tested.tsv")))
    b <- readRich(p)
    expect_round_trip(b, legacy)
    expect_warning(tt <- testedTerms(b), "tested table unavailable")
    expect_null(tt)
    expect_identical(enrichmentInfo(b)$schema_version, 1L)
    # Read back, the object gains the slot prototype; the warning must not
    # then claim it was built by as(), richDAVID() or new().
    msg <- tryCatch(testedTerms(b), warning = conditionMessage)
    expect_match(msg, "no tested table was recorded for this object")
    expect_no_match(msg, "as\\(\\)|richDAVID|new\\(\\)|outside an analysis")
  }
})

test_that("readRich restores a GSEA result that reported nothing", {
  # Catches an empty GSEA result rejected or reloaded with retyped columns.
  skip_if_not_installed("fgsea")
  g <- .fgsea_run(minSize = 3, maxSize = 9, pvalue = 1e-12, padj = NULL)
  expect_identical(nrow(g@result), 0L)
  expect_gt(nrow(testedTerms(g)), 0)
  expect_round_trip(round_trip(g), g)
})

test_that("a filtered view round-trips with the original analysis record", {
  # Catches a derived view stored with a pruned tested table or rewritten
  # metadata (AC-005 I/O leg).
  r <- .fora_reference()
  f <- dplyr::filter(r, Annot != "T1")
  expect_lt(nrow(f@result), nrow(r@result))
  b <- round_trip(f)
  expect_round_trip(b, f)
  expect_identical(testedTerms(b), testedTerms(r))
  expect_identical(enrichmentInfo(b), enrichmentInfo(r))
})

test_that("overwrite = TRUE replaces a previous writeRich folder", {
  # Catches an overwrite that leaves the earlier result's files behind (a stale
  # tested.tsv would make the new folder unreadable).
  w <- io_dir()
  writeRich(.fora_reference(), w)
  x <- methods::as(.fora_reference(pvalue = 1)@result, "richResult")
  expect_no_error(writeRich(x, w, overwrite = TRUE))
  expect_false(file.exists(file.path(w, "tested.tsv")))
  expect_round_trip(readRich(w), x)
})

# ---- readRich: refusals (AC-013) ---------------------------------------------

test_that("readRich refuses a folder it cannot read exactly", {
  # Catches a newer, older or foreign folder read as if it were understood, a
  # truncated or hand-edited table returned as a valid result, and a corrupt
  # number silently read as NA.
  x <- .fora_reference(pvalue = 1)
  expect_identical(nrow(x@result), 2L)
  base <- io_dir()
  writeRich(x, base)
  # The edits below, not the rewriting itself, are what readRich refuses.
  same <- copy_folder(base)
  edit_manifest(same, identity)
  expect_round_trip(readRich(same), x)
  changed <- function(f) {
    p <- copy_folder(base)
    f(p)
    p
  }
  set_field <- function(field, value) {
    changed(function(p) edit_manifest(p, function(m) {
      m[field] <- list(value)
      if (is.null(value)) m[[field]] <- NULL
      m
    }))
  }
  set_schema <- function(value) {
    changed(function(p) edit_metadata(p, function(meta) {
      meta["schema_version"] <- list(value)
      if (is.null(value)) meta$schema_version <- NULL
      meta
    }))
  }

  expect_error(readRich(set_field("format_version", 3L)),
               "format version 3 .*newer.*supported: 1, 2\\)")
  expect_error(readRich(set_field("format_version", 0L)),
               "format version 0 .*older.*supported: 1, 2\\)")
  expect_error(readRich(set_field("format_version", "1")),
               "format version \"1\" .*not a whole number.*supported: 1, 2\\)")
  expect_error(readRich(set_field("format_version", 1.5)),
               "format version 1.5 .*not a whole number.*supported: 1, 2\\)")
  expect_error(readRich(set_field("format_version", NULL)),
               "no format version .*supported: 1, 2\\)")
  expect_error(readRich(set_schema(3L)),
               "schema_version 3 .*newer.*supported: 1, 2\\)")
  expect_error(readRich(set_schema("2")),
               "schema_version \"2\" .*not a whole number.*supported: 1, 2\\)")
  expect_error(readRich(set_schema(NULL)),
               "no metadata schema_version .*supported: 1, 2\\)")
  expect_error(readRich(set_field("format", "other")),
               "format \"other\".*richR-result")
  expect_error(readRich(set_field("class", "Annot")),
               "class \"Annot\".*richResult, GSEAResult")

  expect_error(readRich(changed(function(p) unlink(file.path(p, "tested.tsv")))),
               "tested.tsv.*missing")
  untested <- io_dir()
  writeRich(methods::as(x@result, "richResult"), untested)
  extra <- copy_folder(untested)
  file.copy(file.path(base, "tested.tsv"), extra)
  expect_error(readRich(extra), "tested.tsv.*not listed")
  extra <- copy_folder(untested)
  file.copy(file.path(base, "detail.tsv"), extra)
  expect_error(readRich(extra), "detail.tsv.*not listed")

  expect_error(readRich(changed(function(p) {
    edit_lines(p, "result.tsv", function(l) l[-length(l)])
  })), "result.tsv.* 1 data rows.*records 2")
  expect_error(readRich(changed(function(p) {
    edit_lines(p, "result.tsv", function(l) {
      l[1] <- sub("Pvalue", "Pval", l[1], fixed = TRUE)
      l
    })
  })), "header of result.tsv.*'Pvalue'")
  set_cell <- function(row, column, value) {
    changed(function(p) edit_lines(p, "result.tsv", function(l) {
      cells <- strsplit(l[row + 1L], "\t", fixed = TRUE)[[1]]
      cells[match(column, names(x@result))] <- value
      l[row + 1L] <- paste(cells, collapse = "\t")
      l
    }))
  }
  expect_error(readRich(set_cell(2L, "Pvalue", "0.05x")),
               "result.tsv.*column 'Pvalue', row 2")
  expect_error(readRich(set_cell(1L, "Term", "a\\qb")),
               "result.tsv.*column 'Term', row 1.*unknown escape")
  # A double cell is read only in the forms writeRich() writes: text that
  # as.numeric() would also accept (spaces, hex, other spellings) is refused.
  for (cell in c(" 0x1p-3 ", "0x10", "1e5 ", "+1", "1E5", ".5", "Infinity",
                 "nan", "inf", "1e+999", "-1e+999", "1e-999", "-2.5e-400",
                 paste0("0.", strrep("0", 400), "1"))) {
    expect_error(readRich(set_cell(1L, "Pvalue", cell)),
                 "result.tsv.*column 'Pvalue', row 1.*does not parse as double",
                 info = cell)
  }
  # A data row with a field too many or too few is named by its own row.
  nc <- ncol(x@result)
  edit_row <- function(row, f) {
    changed(function(p) edit_lines(p, "result.tsv", function(l) {
      l[row + 1L] <- f(l[row + 1L])
      l
    }))
  }
  expect_error(readRich(edit_row(2L, function(l) paste0(l, "\textra"))),
               sprintf("result.tsv in '.*', row 2: %d fields found, %d expected",
                       nc + 1L, nc))
  expect_error(readRich(edit_row(1L, function(l) sub("\t[^\t]*$", "", l))),
               sprintf("result.tsv in '.*', row 1: %d fields found, %d expected",
                       nc - 1L, nc))
  expect_error(readRich(changed(function(p) unlink(file.path(p, "manifest.json")))),
               "manifest.json")
})

test_that("readRich decodes only plain data from the manifest", {
  # Catches readRich() rebuilding functions, calls, environments, namespaces,
  # serialized bytes, S4 objects or classed values from a folder it did not
  # write (jsonlite::unserializeJSON() rebuilds all of them, and loads the
  # package a namespace value names, while the folder is read).
  x <- .fora_reference(pvalue = 1)
  base <- io_dir()
  writeRich(x, base)
  node_of <- function(value) {
    jsonlite::fromJSON(jsonlite::serializeJSON(value), simplifyVector = FALSE)
  }
  # A namespace not yet loaded; serializing it would load it, so its node is
  # written by hand.
  spare <- setdiff(c("splines", "stats4", "tcltk"), loadedNamespaces())[1L]
  ns_node <- list(type = "namespace",
                  attributes = structure(list(), names = character(0)),
                  value = list(name = if (is.na(spare)) "splines" else spare))
  inject <- function(node, where = "metadata") {
    p <- copy_folder(base)
    edit_manifest(p, function(m) {
      tree <- m[[where]]
      labels <- vapply(tree$attributes$names$value, as.character, character(1))
      k <- if (where == "metadata") length(labels) + 1L else match("organism", labels)
      tree$value[[k]] <- node
      tree$attributes$names$value[[k]] <- if (where == "metadata") "injected" else "organism"
      m[[where]] <- tree
      m
    })
    p
  }
  before <- loadedNamespaces()
  cases <- list(closure = node_of(function() system("echo PWNED")),
                language = node_of(quote(system("echo PWNED"))),
                environment = node_of(new.env()),
                namespace = ns_node,
                raw = node_of(as.raw(1:3)),
                builtin = node_of(sum),
                S4 = node_of(methods::new("Annot")))
  for (type in names(cases)) {
    expect_error(readRich(inject(cases[[type]])),
                 sprintf("manifest.json of .*type '%s' at metadata\\$injected", type),
                 info = type)
  }
  expect_error(readRich(inject(node_of(as.Date("2026-01-01")))),
               "manifest.json of .*attribute 'class' at metadata\\$injected")
  expect_error(readRich(inject(node_of(function() 1), "slots")),
               "manifest.json of .*type 'closure' at slots\\$organism")
  expect_identical(setdiff(loadedNamespaces(), before), character(0))
  # The plain manifest the edits started from still reads exactly.
  expect_round_trip(readRich(base), x)
})

# SPEC-RICHR-009: the collection itself must not bypass the typed-data guard.
test_that("manifest attributes need unique names and every value is checked", {
  x <- .fora_reference()
  base <- io_dir(); writeRich(x, base)
  nm <- read_manifest(base)$metadata$attributes$names
  for (attrs in list(unname(list(nm)), stats::setNames(list(nm), ""),
                    stats::setNames(list(nm, nm), c("names", "names")))) {
    p <- copy_folder(base)
    edit_manifest(p, function(m) { m$metadata$attributes <- attrs; m })
    expect_error(readRich(p), "manifest.json of .*plain data")
  }
  expect_round_trip(readRich(base), x)
})

test_that("namespace nodes hidden in unnamed or duplicate attributes never load a package", {
  # Fresh child processes keep this assertion meaningful regardless of the
  # namespaces loaded by earlier tests in the suite.
  x <- .fora_reference()
  base <- io_dir(); writeRich(x, base)
  script <- tempfile(fileext = ".R")
  writeLines(c(
    paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
    'library(richR)',
    'args <- commandArgs(TRUE); p <- args[1]; mode <- args[2]',
    'ns <- setdiff(c("splines", "tcltk", "stats4", "compiler"), loadedNamespaces())[1]',
    'stopifnot(!is.na(ns))',
    'm <- jsonlite::fromJSON(file.path(p,"manifest.json"), simplifyVector=FALSE)',
    'evil <- list(type="namespace",attributes=structure(list(),names=character()),value=list(name=ns,version="0"))',
    'done <- FALSE',
    'inject <- function(n) {',
    ' if (done || !is.list(n)) return(n)',
    ' if (identical(n$type,"character")) {',
    '  plain <- list(type="character",attributes=structure(list(),names=character()),value=list("x"))',
    '  n$attributes <- if (mode=="unnamed") list(evil) else setNames(list(plain,evil),c("names","names"))',
    '  done <<- TRUE; return(n)',
    ' }',
    ' if (is.list(n$value)) n$value <- lapply(n$value,inject)',
    ' n',
    '}',
    'm$metadata <- inject(m$metadata); stopifnot(done)',
    'writeLines(jsonlite::toJSON(m,auto_unbox=TRUE,null="null",digits=I(17)),file.path(p,"manifest.json"))',
    'e <- tryCatch({readRich(p); "NO ERROR"},error=conditionMessage)',
    'stopifnot(grepl("manifest.json of .*plain data",e), !ns %in% loadedNamespaces())',
    'cat("NAMESPACE_GUARD_OK\\n")'
  ), script)
  for (mode in c("unnamed", "duplicate")) {
    p <- copy_folder(base)
    out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script), shQuote(p), mode), stdout = TRUE, stderr = TRUE))
    expect_null(attr(out, "status"), info = paste(mode, paste(out, collapse = "\n")))
    expect_true("NAMESPACE_GUARD_OK" %in% out, info = mode)
  }
})

test_that("manifest text is JSON itself, never a file or URL to follow", {
  x <- .fora_reference()
  a <- io_dir(); b <- io_dir()
  writeRich(x, a); writeRich(x, b)
  for (text in c(normalizePath(file.path(b, "manifest.json")), "http://127.0.0.1:9/private")) {
    writeLines(text, file.path(a, "manifest.json"))
    expect_warning(expect_error(readRich(a), "manifest.json"), NA)
  }
  writeLines(normalizePath(file.path(b, "manifest.json")), file.path(a, "manifest.json"))
  writeLines("user data", file.path(a, "result.tsv"))
  before <- dir_state(a)
  expect_error(writeRich(x, a, overwrite = TRUE), "not written by writeRich")
  expect_identical(dir_state(a), before)
})

test_that("extreme doubles retain their bits including negative zero", {
  x <- fs22()
  vals <- c(5e-324, -5e-324, .Machine$double.xmax, -.Machine$double.xmax,
            .Machine$double.xmin, 0, -0, Inf, -Inf, NaN, NA, 1e-300,
            1 / 3, 1 / 7, pi, 0.1, 1 - .Machine$double.eps)
  x@result$Pvalue <- rep_len(vals, nrow(x@result))
  y <- round_trip(x)
  expect_identical(y@result$Pvalue, x@result$Pvalue)
  expect_identical(writeBin(y@result$Pvalue, raw()), writeBin(x@result$Pvalue, raw()))
  expect_identical(1 / y@result$Pvalue[7], -Inf)
})

test_that("overwrite treats wildcard characters literally and removes stale tested tables", {
  x <- .fora_reference()
  bare <- methods::as(x@result, "richResult")
  base <- io_dir(); dir.create(base)
  pairs <- list(c("out[1]", "out1"))
  # Windows cannot create '*'; the bracket case still exercises literal paths.
  if (.Platform$OS.type != "windows") pairs <- c(pairs, list(c("star*", "starX")))
  for (pair in pairs) {
    target <- file.path(base, pair[1]); sibling <- file.path(base, pair[2])
    dir.create(sibling)
    for (f in c("manifest.json", "result.tsv", "tested.tsv"))
      writeLines("precious", file.path(sibling, f))
    before <- dir_state(sibling)
    writeRich(x, target)
    writeRich(bare, target, overwrite = TRUE)
    expect_identical(dir_state(sibling), before)
    expect_false(file.exists(file.path(target, "tested.tsv")))
    expect_round_trip(readRich(target), bare)
  }
})

test_that("overwrite checks folder writability before deleting the old result", {
  x <- .fora_reference()
  p <- io_dir(); writeRich(x, p)
  before <- dir_state(p)
  Sys.chmod(p, "0555")
  on.exit(Sys.chmod(p, "0755"), add = TRUE)
  if (file.access(p, 2) == 0L) skip("folder remains writable on this filesystem/user")
  expect_error(writeRich(methods::as(x@result, "richResult"), p, overwrite = TRUE),
               "writeRich.*(writable|remove)")
  expect_identical(dir_state(p), before)
  expect_round_trip(readRich(p), x)
})

test_that("a list cell must itself be JSON even when a matching local file exists", {
  d <- data.frame(pathway = "P", Term = "term", Pvalue = 0.1)
  d$leadingEdge <- list("safe")
  x <- methods::new("GSEAResult", result = d, input = c(gene = 1))
  p <- io_dir(); writeRich(x, p)
  edit_lines(p, "result.tsv", function(l) {
    cells <- strsplit(l[2], "\t", fixed = TRUE)[[1]]
    cells[match("leadingEdge", names(d))] <- "[x"
    l[2] <- paste(cells, collapse = "\t")
    l
  })
  cwd <- io_dir(); dir.create(cwd)
  old <- setwd(cwd); on.exit(setwd(old), add = TRUE)
  writeLines('["LEAK"]', "[x")
  expect_error(readRich(p), "result.tsv.*column 'leadingEdge', row 1.*does not parse as list")
})
test_that("overwrite stops when a known entry cannot be removed", {
  x <- .fora_reference()
  p <- io_dir(); writeRich(x, p)
  on.exit(unlink(p, recursive = TRUE, expand = FALSE), add = TRUE)
  manifest <- file.path(p, "manifest.json")
  old_manifest <- readBin(manifest, "raw", n = file.info(manifest)$size)
  replacement <- methods::as(x@result, "richResult")
  unlink(file.path(p, "detail.tsv"))
  dir.create(file.path(p, "detail.tsv"))
  marker <- file.path(p, "detail.tsv", "keep.txt")
  writeLines("keep", marker)
  expect_error(writeRich(replacement, p, overwrite = TRUE), "writeRich.*could not remove")
  expect_identical(readLines(marker), "keep")
  # The original manifest may remain or be removed, but cannot be replaced.
  if (file.exists(manifest)) {
    expect_identical(readBin(manifest, "raw", n = file.info(manifest)$size), old_manifest)
  }
})
test_that("overwrite expands a home-relative path without expanding its wildcards", {
  home <- path.expand("~")
  if (file.access(home, 2L) != 0L) skip("home directory is not writable")
  root <- tempfile("rich-io-home-", tmpdir = home)
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, expand = FALSE), add = TRUE)
  p <- paste0("~/", basename(root), "/out[1]")
  sibling <- file.path(root, "out1")
  dir.create(sibling)
  keep <- file.path(sibling, "tested.tsv")
  writeLines("precious", keep)
  x <- .fora_reference()
  writeRich(x, p)
  expect_true(file.exists(file.path(path.expand(p), "tested.tsv")))
  bare <- methods::as(x@result, "richResult")
  expect_identical(withVisible(writeRich(bare, p, overwrite = TRUE)),
                   list(value = p, visible = FALSE))
  expect_identical(readLines(keep), "precious")
  expect_false(file.exists(file.path(path.expand(p), "tested.tsv")))
  expect_identical(readRich(p), bare)
})

test_that('readback rejects a self-consistent table without its required identifier', {
  p <- io_dir(); writeRich(.fora_reference(),p)
  edit_lines(p,'result.tsv',function(l) {l[1] <- sub('^Annot','unknown',l[1]);l})
  edit_manifest(p,function(m) {m$tables[[1]]$columns[[1]]$name <- 'unknown';m})
  expect_error(readRich(p),'result.*lacks.*Annot|result.*required.*Annot')
})

test_that('legacy format one is explicit about unavailable checksums', {
  r <- .fora_reference(); p <- io_dir(); writeRich(r,p)
  edit_manifest(p,function(m) {
    m$format_version <- 1L
    m$tables <- lapply(m$tables,function(t) {t$sha256 <- NULL;t})
    m
  })
  expect_warning(z <- readRich(p),'checksum.*unavailable|unavailable.*checksum')
  expect_round_trip(z,r)
})
