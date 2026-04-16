# richR Updates

## Version 0.1.3 (2026-04-16)

### Bug Fixes

- **`richUpset()` bar-to-dot alignment**: The intersection bar chart
  and the dot matrix previously used a nested `plot_grid()` layout
  where each row laid out the left gutter and the main plot
  independently. Differences in y-axis label widths (bar-chart tick
  labels vs. set-name labels) produced distinct panel start positions,
  causing the bars to drift out of alignment with the dots beneath
  them. Panels are now pre-aligned with `cowplot::align_plots()` and
  placed at explicit coordinates with `ggdraw()` + `draw_plot()` so
  the bar chart and matrix always share identical x-ranges regardless
  of label width. The `set.size.show = FALSE` branch additionally
  gained `align = "v", axis = "lr"` to `plot_grid()` for the same
  reason.

### Enhancements

- **`richUpset()` multi-set bar colors**: Previously, bars representing
  multi-set intersections inherited the color of the largest
  participating set. Bars for different intersections could therefore
  share a color, making it hard to tell visually which sets were
  combined. Multi-set bars now use an RGB-blend of the participating
  set colors (computed by the new internal helper `.blend_colors()`).
  Single-set bars still use the set color, and setting
  `main.bar.color` still overrides every bar to a uniform color.

### Internal

- New internal helper `.blend_colors()` in `R/ggupset.R` that averages
  RGB channels of a vector of color specifications.
- NAMESPACE: added `importFrom(cowplot, align_plots)`,
  `importFrom(cowplot, ggdraw)`, `importFrom(cowplot, draw_plot)`,
  `importFrom(grDevices, col2rgb)`, `importFrom(grDevices, rgb)`.

---

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
- New functions (`readGMT`, `buildAnnotFromList`, `getGenes`, `theme_richR`,
  `batchEnrich`) are explicitly exported.
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
improvements, and major refactoring to improve maintainability.
