# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

richR is an R package for functional enrichment analysis (GO, KEGG, Reactome, MSigDB, GSEA) with visualization. Originally developed by Kai Guo (guokai8/richR), now maintained at hurlab/richR.

## Repository Setup

- **origin**: hurlab/richR (main branch) - the maintained version
- **upstream**: guokai8/richR (master branch) - the original

## Build & Check Commands

```bash
# R CMD check (standard R package validation)
R CMD build . && R CMD check richR_*.tar.gz

# Or from within R:
# devtools::check()
# devtools::document()  # regenerate man pages from roxygen2
# devtools::test()      # run testthat tests
```

## Architecture

### S4 Class System

The package uses S4 classes defined in `R/00-AllClasses.R`:

- **`Annot`**: Annotation database (species, anntype, keytype, annot data.frame)
- **`richResult`**: ORA enrichment result (result, detail, gene, cutoffs, metadata)
- **`GSEAResult`**: GSEA result (result, input fold changes, gene, metadata)

S4 generics are in `R/AllGenerics.R`. S4 methods for accessors (head, tail, dim, $, [, result, detail, summary) are in `R/misc.R`.

### Core Analysis Pipeline

1. **Annotation building** (`R/makeAnnot.R`): `buildAnnot()` → species-specific Bioconductor DB → `Annot` object. `buildMSIGDB()` for MSigDB.
2. **ORA** (`R/enrich.R`, `R/richGO.R`, `R/richKEGG.R`): hypergeometric test via C++ (`src/hyper.cpp`). Each has `*_internal()` function + S4 method dispatchers for data.frame and Annot signatures.
3. **GSEA** (`R/richGSEA.R`): wraps `fgsea::fgseaMultilevel()`. Includes `plotGSEA()`, `getPathways()`, `searchPathways()`, `filterPathways()`, `batchGSEAplot()`, `parGSEA()`.
4. **Clustering** (`R/kappa.R`): kappa-score similarity → `.merge_term()` → enrichment scoring.
5. **Visualization**: `R/ggbar.R`, `R/ggdot.R`, `R/ggrich.R` (network), `R/ggnetwork.R`, `R/ggnetmap.R`, `R/ggdot_cluster.R` + `R/ggdot_cluster_gsea.R`, `R/compareResult.R`, `R/ggheatmap.R`.

### C++ (Rcpp) Components (`src/`)

- `hyper.cpp`: Vectorized hypergeometric p-value computation
- `sf.cpp`: Fast split-factor (named list from two-column data.frame)
- `name_table.cpp`: Fast tabulation of named list lengths
- `reverse_list.cpp`: C++ implementation of key-value reversal
- `unique.cpp`: Fast unique element extraction

### dplyr Integration

`R/dplyr-*.R` files define S3 methods for dplyr verbs on richResult/GSEAResult/Annot. Re-exports in `R/dplyr_exportR.R`.

### Key Data Files (`data/`)

- `godata.rdata`: Prebuilt GO offspring data
- `kegg.rdata`: Built-in KEGG pathway annotations
- `path.rda`: KEGG pathway hierarchy (Level1/Level2/Level3)
- `module.rda`: KEGG module annotations
- `pathway.rda`: Pathway lookup data

## Version Policy

- **Minor version** (0.X.0): Feature additions, significant doc rewrites, bug fix batches
- **Patch version** (0.0.X): Single feature additions, small fixes
- All version changes must be logged in `UPDATES.md`

## Important Conventions

- NAMESPACE is roxygen2-generated: edit roxygen comments in R files, then run `devtools::document()`
- The `exportPattern("^[[:alpha:]]+")` in NAMESPACE exports all public functions; individual `@export` tags are also present
- Species lookup uses common names ("human", "mouse", etc.) mapped to Bioconductor package names in `.getdb()` and KEGG codes in `.getspeices()`
- The package sets `options(stringsAsFactors = FALSE)` on load (`R/zzz.R`)
