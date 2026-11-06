# mapGenes() (SPEC-RICHR-005 REQ-001 .. REQ-005): explicit identifier
# conversion with one recorded decision per input entry. Fixture F-MAP is
# inline; only the OrgDb leg needs org.Hs.eg.db (a Suggests package).

fmap_table <- function() {
  data.frame(
    from = c("a", "a", "b", "b", "d", "e", "g", "h", "k", paste0("z", 1:10)),
    to = c("A1", "A1", "B1", "B2", "A1", "E1", "E1", "H1", NA, paste0("Z", 1:10)),
    stringsAsFactors = FALSE
  )
}
fmap_query <- c("a", "b", "c", "d", NA, "", "e", "a")
fmap_ranking <- function() {
  s <- c(2, 2, 1.5, -1, 5, 0.7, 1, NaN, 4, 0.7, NaN, 3, seq(0.1, 1.0, by = 0.1))
  names(s) <- c("a", "d", "b", "e", "g", "h", "c", "f", NA, "h", "h", "k", paste0("z", 1:10))
  s
}
mapping_record <- function(x) attr(x, "richR_mapping", exact = TRUE)

# The canonical text of design.md § 5, built here independently of the package.
idmap_canonical <- function(from, to) {
  esc <- function(v) {
    v <- enc2utf8(as.character(v))
    v <- gsub("\\", "\\\\", v, fixed = TRUE)
    v <- gsub("\t", "\\t", v, fixed = TRUE)
    v <- gsub("\n", "\\n", v, fixed = TRUE)
    v <- gsub("\r", "\\r", v, fixed = TRUE)
    v[is.na(v)] <- "\\N"
    v
  }
  paste0("richR-idmap-v1\n", paste0(esc(from), "\t", esc(to), "\n", collapse = ""))
}

# Any attribute other than names makes M4's readRich() refuse the record.
plain_attributes_only <- function(x) {
  extra <- setdiff(names(attributes(x)), "names")
  if (length(extra)) return(FALSE)
  if (is.list(x)) return(all(vapply(x, plain_attributes_only, logical(1))))
  TRUE
}

test_that("AC-001: one-to-many ids are dropped visibly and every entry keeps one decision", {
  # Catches: a one-to-many identifier converted silently (idconvert(), L-M1),
  # a duplicated source row read as ambiguous, missing or blank entries dropped
  # without a row, and an output whose length disagrees with its record.
  tab <- fmap_table()
  v <- mapGenes(fmap_query, tab, "from", "to")
  expect_identical(as.vector(v), c("A1", "E1"))

  rec <- mapping_record(v)
  expect_true(plain_attributes_only(rec))
  expect_identical(names(rec$decisions),
                   c("position", "input_id", "candidates", "target", "decision", "reason"))
  expect_null(attr(rec$decisions, "class"))
  expect_null(attr(rec$decisions, "row.names"))

  d <- mappingAudit(v)
  expect_s3_class(d, "data.frame")
  expect_identical(names(d), names(rec$decisions))
  expect_true(is.list(d$candidates))
  expect_identical(.row_names_info(d), -8L)
  expect_identical(d$position, 1:8)
  expect_identical(d$input_id, fmap_query)
  expect_identical(d$decision, c("mapped", "ambiguous", "unmapped", "merged",
                                 "invalid", "invalid", "mapped", "merged"))
  expect_identical(d$target, c("A1", NA, NA, "A1", NA, NA, "E1", "A1"))
  expect_identical(d$candidates[[1]], "A1")
  expect_identical(d$candidates[[2]], c("B1", "B2"))
  expect_identical(d$candidates[[3]], character(0))
  expect_identical(length(unique(d$target[d$decision == "mapped"])), length(v))

  first <- mapGenes(fmap_query, tab, "from", "to", multiple = "first")
  expect_identical(as.vector(first), c("A1", "B1", "E1"))
  df <- mappingAudit(first)
  expect_identical(df$decision[[2]], "mapped")
  expect_match(df$reason[[2]], "2 candidates")

  # Under "drop" the result does not depend on the source's row order; under
  # "first" the row order decides, and visibly so.
  rev_tab <- tab[rev(seq_len(nrow(tab))), ]
  v_rev <- mapGenes(fmap_query, rev_tab, "from", "to")
  expect_identical(as.vector(v_rev), as.vector(v))
  expect_identical(mappingAudit(v_rev)$decision, d$decision)
  expect_identical(mappingAudit(v_rev)$target, d$target)
  first_rev <- mapGenes(fmap_query, rev_tab, "from", "to", multiple = "first")
  expect_identical(mappingAudit(first_rev)$target[[2]], "B2")
})

test_that("AC-002: a ranking keeps scores, merges equal collisions and never picks by score", {
  # Catches: duplicate names richGSEA() rejects (L-M8), a collision resolved by
  # keeping the stronger score, an invalid entry poisoning a valid target, a
  # missing target treated as a candidate, a non-finite score passed on.
  s <- fmap_ranking()
  r <- mapGenes(s, fmap_table(), "from", "to")
  want <- c(A1 = 2, H1 = 0.7, stats::setNames(seq(0.1, 1.0, by = 0.1), paste0("Z", 1:10)))
  expect_identical(as.vector(r), unname(want))
  expect_identical(names(r), names(want))

  rec <- mapping_record(r)
  expect_true(plain_attributes_only(rec))
  expect_identical(names(rec$decisions),
                   c("position", "input_id", "candidates", "target", "decision",
                     "reason", "score"))
  expect_identical(rec$decisions$score, unname(as.double(s)))

  d <- mappingAudit(r)
  expect_identical(d$decision[1:12],
                   c("mapped", "merged", "ambiguous", "conflict", "conflict", "mapped",
                     "unmapped", "invalid", "invalid", "merged", "invalid", "unmapped"))
  expect_true(all(d$decision[13:22] == "mapped"))
  expect_false("E1" %in% names(r))
  expect_true(all(is.na(d$target[d$decision %in% c("conflict", "invalid", "unmapped", "ambiguous")])))

  skip_if_not_installed("fgsea")
  ann <- data.frame(GeneID = names(r), Term = rep(c("T1", "T2"), length.out = length(r)))
  expect_no_error(suppressWarnings(
    richGSEA(r, ann, minSize = 2, maxSize = 50, padj = NULL, pvalue = 1, seed = 1)))
})

test_that("AC-003: the record names its source and a reproducible table fingerprint", {
  # Catches: a conversion whose source cannot be recovered, and a fingerprint
  # not reproducible from its stated canonical form or blind to a changed or
  # reordered mapping.
  tab <- fmap_table()
  v <- mapGenes(fmap_query, tab, "from", "to")
  src <- mapping_record(v)$source
  expect_identical(src$kind, "table")
  expect_identical(src$from, "from")
  expect_identical(src$to, "to")
  expect_identical(src$multiple, "drop")
  expect_identical(src$n_rows, nrow(tab))
  expect_identical(src$fingerprint$algorithm, "sha256")
  expect_identical(src$fingerprint$canonicalization, "richR-idmap-v1")
  expect_identical(src$fingerprint$hash,
                   digest::digest(idmap_canonical(tab$from, tab$to), algo = "sha256",
                                  serialize = FALSE))
  expect_identical(mapGenes(fmap_query, tab, "from", "to"), v)

  odd <- data.frame(from = c("x\ty", "p\\q", "m"), to = c("T1", "T2", NA))
  fo <- mapping_record(mapGenes(c("x\ty", "m"), odd, "from", "to"))$source$fingerprint
  expect_identical(fo$hash, digest::digest(idmap_canonical(odd$from, odd$to),
                                           algo = "sha256", serialize = FALSE))

  hash_of <- function(t) mapping_record(mapGenes("a", t, "from", "to"))$source$fingerprint$hash
  changed <- tab
  changed$to[6] <- "E2"
  swapped <- tab[c(2, 1, 4, 3, 5:nrow(tab)), ]
  expect_false(identical(hash_of(changed), src$fingerprint$hash))
  expect_false(identical(hash_of(swapped), src$fingerprint$hash))
  # A factor column is read as its labels: same decisions, same fingerprint.
  ftab <- data.frame(from = factor(tab$from), to = factor(tab$to))
  expect_identical(hash_of(ftab), src$fingerprint$hash)

  skip_if_not_installed("org.Hs.eg.db")
  db <- org.Hs.eg.db::org.Hs.eg.db
  keys <- c("23", "1", "999999999")
  o <- mapGenes(keys, db, "ENTREZID", "ENSEMBL")
  od <- mappingAudit(o)
  all23 <- suppressMessages(AnnotationDbi::mapIds(db, keys = "23", column = "ENSEMBL",
                                                  keytype = "ENTREZID", multiVals = "list"))
  expect_identical(od$decision, c("ambiguous", "mapped", "unmapped"))
  expect_identical(od$candidates[[1]], unname(all23[["23"]]))
  of <- mapGenes(keys, db, "ENTREZID", "ENSEMBL", multiple = "first")
  expect_identical(mappingAudit(of)$target[[1]],
                   unname(suppressMessages(idconvert("human", "23", "ENTREZID", "ENSEMBL"))))
  osrc <- mapping_record(o)$source
  expect_true(plain_attributes_only(mapping_record(o)))
  expect_identical(osrc$kind, "OrgDb")
  expect_identical(osrc$package, "org.Hs.eg.db")
  expect_identical(osrc$version, as.character(utils::packageVersion("org.Hs.eg.db")))
  expect_identical(names(osrc$metadata), c("name", "value"))
  expect_identical(as.data.frame(osrc$metadata, stringsAsFactors = FALSE),
                   AnnotationDbi::metadata(db))
  expect_identical(mapGenes(keys, "human", "ENTREZID", "ENSEMBL"), o)
  expect_identical(mapGenes(keys, db, "ENTREZID", "ENSEMBL"), o)
})

test_that("AC-004: malformed calls error naming the argument; an all-dropped call warns", {
  # Catches: a malformed call turned into an empty or wrong mapping, a double
  # identifier column matched as "1e+05" (L-N1), an absent OrgDb package
  # surfacing obscurely, and an all-dropped conversion hiding its decisions.
  tab <- fmap_table()
  expect_error(mapGenes(list("a"), tab, "from", "to"), "`x`.*list")
  expect_error(mapGenes(character(0), tab, "from", "to"), "`x`.*empty")
  expect_error(mapGenes(c(1, 2), tab, "from", "to"), "`x`.*without names")
  expect_error(mapGenes("a", tab, "from", "from"), "`from`.*`to`.*\"from\"")
  expect_error(mapGenes("a", tab, "from", "missing_col"), "`to`.*missing_col")
  expect_error(mapGenes("a", tab, "from", "to", multiple = "all"), "`multiple`.*all")
  expect_error(mapGenes("100000", data.frame(from = c(100000, 23), to = c("A", "B")),
                        "from", "to"), "`from`.*numeric")
  expect_error(mapGenes("a", "no_such_species", "ENTREZID", "ENSEMBL"),
               "`mapping`.*no_such_species")
  local_mocked_bindings(.getdb = function(...) "richRnotInstalledPkg.db", .package = "richR")
  expect_error(mapGenes("1", "human", "ENTREZID", "ENSEMBL"),
               "`mapping`.*richRnotInstalledPkg\\.db")

  expect_warning(none <- mapGenes(c("zz", NA, "b"), tab, "from", "to"),
                 "invalid: 1.*unmapped: 1.*ambiguous: 1")
  expect_identical(as.vector(none), character(0))
  expect_identical(mappingAudit(none)$decision, c("unmapped", "invalid", "ambiguous"))
})
