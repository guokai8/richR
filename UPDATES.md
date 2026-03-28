# richR Updates

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
