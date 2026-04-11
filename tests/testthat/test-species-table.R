test_that(".species_table returns valid data.frame", {
  tbl <- richR:::.species_table()
  expect_s3_class(tbl, "data.frame")
  expect_true(all(c("species", "dbname", "kegg", "msigdb", "reactome") %in% colnames(tbl)))
  expect_true(nrow(tbl) >= 20)
  expect_false(any(duplicated(tbl$species)))
  expect_false(any(is.na(tbl$species)))
  expect_false(any(is.na(tbl$dbname)))
})

test_that(".getdb returns correct OrgDb for known species", {
  expect_equal(richR:::.getdb("human"), "org.Hs.eg.db")
  expect_equal(richR:::.getdb("mouse"), "org.Mm.eg.db")
  expect_equal(richR:::.getdb("zebrafish"), "org.Dr.eg.db")
  expect_equal(richR:::.getdb("yeast"), "org.Sc.sgd.db")
  expect_null(richR:::.getdb("unicorn"))
})

test_that(".getspeices returns correct KEGG codes", {
  expect_equal(richR:::.getspeices("human"), "hsa")
  expect_equal(richR:::.getspeices("mouse"), "mmu")
  expect_equal(richR:::.getspeices("fly"), "dme")
  expect_null(richR:::.getspeices("unicorn"))
})

test_that(".getmsig returns correct scientific names", {
  expect_equal(richR:::.getmsig("human"), "Homo sapiens")
  expect_equal(richR:::.getmsig("mouse"), "Mus musculus")
  expect_equal(richR:::.getmsig("fly"), "Drosophila melanogaster")
  expect_null(richR:::.getmsig("unicorn"))
})

test_that(".getrodbname returns correct Reactome species", {
  expect_equal(richR:::.getrodbname("human"), "Homo sapiens")
  expect_equal(richR:::.getrodbname("mouse"), "Mus musculus")
  expect_equal(richR:::.getrodbname("arabidopsis"), "Arabidopsis thaliana")
  expect_null(richR:::.getrodbname("unicorn"))
})

test_that("showData returns consistent data.frame", {
  dd <- showData()
  expect_s3_class(dd, "data.frame")
  expect_true(all(c("dbname", "species") %in% colnames(dd)))
  expect_true("org.Hs.eg.db" %in% dd$dbname)
  expect_true("human" %in% dd$species)
})
