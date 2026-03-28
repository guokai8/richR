# richR Updates

## Version 0.1.1.1 (2026-03-28)

### Bug Fixes

- **UpSet plot bar-to-dot alignment**: Bars and dots were progressively
  misaligned because the bar chart and dot matrix panels had different left
  margins. Fixed by using `align_plots()` + `ggdraw()`/`draw_plot()` for
  pixel-accurate panel alignment, replacing the nested `plot_grid()` layout.
- **UpSet multi-set intersection bar colors**: Bars representing multi-set
  intersections now use a blended color (RGB average of participating sets)
  instead of arbitrarily picking the largest set's color. Single-set bars
  still use the set color. Override all bars with `main.bar.color`.

### Documentation

- README UpSet example now uses a self-contained reproducible dataset
  (`set.seed(123)`, 250-gene universe) matching the example figure.
- Fixed citation year (2025 → 2026) and version.

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
