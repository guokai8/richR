# richR Updates

## Version 0.1.2 (2026-03-28)

### New Features

- **`show()` methods** for richResult, GSEAResult, and Annot classes: typing an
  object at the console now prints a concise summary (organism, ontology, gene
  count, top terms) instead of raw S4 slot dump.
- **Input validation**: all enrichment functions now validate gene input (empty
  lists, NAs, duplicates), check parameter ranges (minGSSize < maxGSSize), and
  report gene-annotation mapping rates. GSEA validates for named numeric input
  and removes NA/Inf values.
- **Empty result handling**: enrichment functions now return informative messages
  when no significant terms are found, instead of cryptic errors.
- **`readGMT()`**: import gene sets from GMT files (the community standard format
  used by MSigDB, Enrichr, GSEA) into Annot objects.
- **`buildAnnotFromList()`**: convert named lists of gene vectors into Annot
  objects for custom enrichment analysis.
- **`getGenes()`**: extract genes from a specific enrichment term, properly
  parsing the comma-separated GeneID field.
- **`theme_richR()`**: consistent publication-ready ggplot2 theme for all richR
  visualizations.
- **`batchEnrich()`**: run enrichment on multiple gene lists in one call, with
  progress reporting and combined output ready for `comparedot()`.
- **Unit test infrastructure**: added testthat tests for S4 class construction,
  show methods, input validation, accessor methods, and edge cases.

### Bug Fixes

- **`data()` anti-pattern**: replaced bare `data(kegg)`, `data(path)`,
  `data(module)` calls inside functions with `data(..., envir)` to avoid
  polluting the global environment.

### Internal

- **NAMESPACE cleanup**: removed overly broad `exportPattern("^[[:alpha:]]+")`.
  Internal functions (`enrich_internal`, `richGSEA_internal`, `vec_to_df`, etc.)
  are no longer exported. New functions (`ggupset`, `readGMT`,
  `buildAnnotFromList`, `getGenes`, `theme_richR`, `batchEnrich`) are explicitly
  exported.
- Added verbose `message()` calls to `buildAnnot()` reporting progress and
  database statistics.

---

## Version 0.1.1 (2026-03-28)

### New Features

- **`ggupset()`**: Color-enhanced UpSet plot for visualizing set intersections
  across multiple gene lists. Built entirely with ggplot2 (no dependency on the
  UpSetR package). Inspired by the colored UpSet implementation from the
  katelynhur/covid19VAE project.
  - Per-set coloring of bars, matrix dots, and set labels via `mycol` parameter
  - Customizable intersection line color via `line.color`
  - Toggle per-set vs. uniform dot coloring via `color.by.set`
  - Order intersections by frequency or degree
  - Accepts named lists of genes, richResult objects, or GSEAResult objects
  - Bar labels colored to match their respective bars

### Documentation

- Added UpSet plot section to README.md with examples and visualization image
- Version bump: 0.1.0 -> 0.1.1

---

## Version 0.1.0 (2026-03-28)

This release includes a comprehensive codebase review, bug fixes, documentation
improvements, and updated package metadata. Based on version 0.0.40 from the
original repository (guokai8/richR).

### Bug Fixes

- **Operator precedence bug** in `enrich_internal()`, `richGO_internal()`, and
  `richKEGG_internal()`: `length(as.vector(resultFis$GeneID)>=1)` was evaluating
  the comparison inside `as.vector()` instead of comparing against `length()`.
  Fixed to `length(as.vector(resultFis$GeneID)) >= 1`.
- **Undefined variable** in `richGSEA` S4 methods: `ontology=ontology` referenced
  an undefined variable in the method signature. Changed to `ontology=""`.
- **Slot access after overwrite** in `ggnetplot` for GSEAResult: `object@sep` was
  accessed after `object` was already overwritten to a data.frame. Fixed extraction
  order.
- **Broken `setAs`** from data.frame to richResult: referenced undefined variables
  (`detail`, `pvalue`, `padj.method`, `input`). Fully rewritten with correct
  variable handling.
- **`summary()` crash** for GSEAResult: `table(...)[[2]]` would error when no
  results had padj < 0.05. Replaced with `sum(..., na.rm = TRUE)`.
- **Function name mismatch**: `getann()` called `.get_kgm_dat()` but the actual
  function name was `.get_kgm.data()`. Fixed the call.
- **Duplicate `.clean.char()`** definition in misc.R. Removed the duplicate.
- **Drosophila typo**: `.getmsig()` returned `"rosophila melanogaster"` (missing
  "D"). Fixed to `"Drosophila melanogaster"`.
- **`showData()` mismatch**: Species names and database package names were not
  aligned. Corrected ordering and fixed "sco" to "streptomyces".

### Deprecation Fixes

- Replaced deprecated ggplot2 `size` parameter with `linewidth` in all
  `geom_segment()` and `geom_curve()` calls across `richGSEA.R`,
  `compareResult.R`, `ggdot_cluster.R`, and `ggdot_cluster_gsea.R`.

### Documentation

- **DESCRIPTION**: Updated title, description, and authors. Added URL and
  BugReports fields. Fixed "for for" double-word typo.
- **README.md**: Complete rewrite with comprehensive examples covering ORA, GSEA,
  visualization, multi-group comparison, term clustering, tidyverse integration,
  and species support table. Includes example plot images.
- **Vignette**: Updated title, fixed language errors ("ensemble" to "Ensembl",
  "mutiple" to "multiple", etc.), improved code formatting.
- **roxygen documentation**: Fixed "bulitin" to "builtin", "bulit" to "built",
  "avaliable" to "available", "databse" to "database", "annotaion" to
  "annotation", "detial" to "detail", "infomation" to "information" across all
  source files and man pages.
- **Removed colons from `@param` tags** (R documentation convention).
- **Fixed duplicate `@export @author`** annotations in `idconvert()`.

### Internal

- Added `Authors@R` field in DESCRIPTION with proper person() entries.
- Package maintained at https://github.com/hurlab/richR.

---

## Version 0.0.40 (2025-12-10)

*Original release by Kai Guo (guokai8/richR). Changes prior to this version
are tracked in the original repository.*
