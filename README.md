# richR <img src="man/figures/logo.png" align="right" height="120" />

[![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![](https://img.shields.io/badge/version-0.1.2-green.svg)](https://github.com/guokai8/richR)
[![R-CMD-check](https://github.com/guokai8/richR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/guokai8/richR/actions/workflows/R-CMD-check.yaml)
![](https://img.shields.io/github/languages/code-size/guokai8/richR)
[![DOI](https://zenodo.org/badge/243827597.svg)](https://zenodo.org/badge/latestdoi/243827597)

## Overview

**richR** is an R package for functional enrichment analysis and visualization.
It supports Over-Representation Analysis (ORA) via hypergeometric testing,
Gene Set Enrichment Analysis (GSEA), and kappa-score-based term clustering.
Built-in annotation builders cover GO, KEGG, Reactome, KEGG Module, and MSigDB,
and you can supply custom gene sets via GMT files or named lists.

### Key features

- **Enrichment**: `richGO()`, `richKEGG()`, `enrich()`, `richGSEA()`, `richDAVID()`
- **Clustering**: `richCluster()` groups terms by kappa similarity
- **20+ visualization functions** including bar, dot, lollipop, heatmap, network,
  volcano, circle bar, gene-dot, gene-heat, gene-bar, UpSet, and more
- **dplyr support**: `filter()`, `select()`, `mutate()`, `arrange()`, `group_by()`,
  `slice()`, `summarise()` work directly on enrichment result objects
- **Batch analysis**: `batchEnrich()` runs enrichment across multiple gene lists
- **Comparison**: `compareResult()` and `richCompareDot()` for multi-group comparisons

## Installation

```r
# Install from GitHub
library(devtools)
install_github("guokai8/richR")
```

Bioconductor annotation packages are needed for building annotations:
```r
BiocManager::install(c("org.Hs.eg.db", "GO.db"))       # human GO
BiocManager::install("reactome.db")                      # Reactome (optional)
```

## Quick Start

### 1. Build annotation data

```r
library(richR)

# Check available species
showData()

# Build GO and KEGG annotations
hsago <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "GO")
hsako <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGG")

# KEGG Module
hsakom <- buildAnnot(species = "human", keytype = "SYMBOL", anntype = "KEGGM")

# MSigDB
hsamgi <- buildMSIGDB(species = "human", keytype = "SYMBOL", anntype = "GO")
```

### 2. Run enrichment analysis

```r
# GO enrichment (ORA)
gene <- sample(unique(hsago$GeneID), 1000)
resgo <- richGO(gene, godata = hsago, ontology = "BP")
head(resgo)

# KEGG enrichment
resko <- richKEGG(gene, kodata = hsako, pvalue = 0.05)

# GSEA (requires a named numeric vector)
genelist <- rnorm(1000)
names(genelist) <- sample(unique(hsako$GeneID), 1000)
res_gsea <- richGSEA(genelist, object = hsako)

# DAVID (online)
res_david <- richDAVID(gene, keytype = "ENTREZID", species = "human")
```

### 3. Custom gene sets

```r
# Import from GMT file (MSigDB, Enrichr, etc.)
annot <- readGMT("h.all.v2023.2.Hs.symbols.gmt", species = "human")
res <- enrich(my_genes, annot)

# From a named list
my_sets <- list(
  "Apoptosis"  = c("TP53", "BAX", "BCL2", "CASP3"),
  "Cell Cycle" = c("CDK1", "CDK2", "CCND1", "RB1")
)
annot <- buildAnnotFromList(my_sets)
res <- enrich(my_genes, annot)
```

### 4. Custom annotation with bioAnno

```r
# library(bioAnno)
# fromKEGG(species = "ath")
# athgo <- buildOwn(dbname = "org.ath.eg.db", anntype = "GO")
# athko <- buildOwn(dbname = "org.ath.eg.db", anntype = "KEGG")
# See https://github.com/guokai8/bioAnno for details
```

## Visualization

richR provides a wide range of plotting functions. All accept `richResult`,
`GSEAResult`, or `data.frame` objects and support saving to file via
`filename`, `width`, and `height` arguments. Each function has a `rich*`
primary name and backward-compatible `gg*` alias.

### Basic plots

```r
# Bar plot and dot plot
richBar(resgo, top = 20, usePadj = FALSE)
richDot(resko, top = 10, usePadj = FALSE)

# Lollipop plot (x-axis = RichFactor, color = significance)
richLollipop(resko, top = 20)

# Circular bar plot
richCircle(resko, top = 15)
```

### Gene-level plots

```r
# Gene-term dot plot (color by fold change or significance)
richGeneDot(resgo, fc = my_foldchanges)
richGeneDot(resgo)                          # no fc: color by -log10(Pvalue)

# Gene-term heatmap
richGeneHeat(resgo, fc = my_foldchanges)
richGeneHeat(resgo, usePadj = FALSE,        # custom colors
             low = "deepskyblue", high = "darkorange", mid = "white")

# Gene bar plot (x-axis = RichFactor, gene labels aligned left)
richGeneBar(resgo, top = 10)
```

### GSEA-specific plots

```r
# NES bar plot
richNES(res_gsea, top = 20)

# GSEA running score curve
richGSEAcurve(res_gsea, term = "hsa04110")

# Volcano plot (NES vs. -log10 pvalue)
richVolcano(res_gsea)

# ECDF plot
richECDF(res_gsea)

# Scatter plot
richScatter(res_gsea)
```

### Network and similarity plots

```r
# Gene-concept network
richNetplot(resko, top = 20)

# Term similarity network
richNetwork(resgo, top = 20, weightcut = 0.01)

# Term similarity heatmap
richTermSim(resko, top = 20)

# Combined network map (multiple enrichment results)
richNetmap(list(resgo, resko), top = 50)
```

### Cluster and comparison plots

```r
# Kappa-based clustering
resc <- richCluster(resgo)
richClusterDot(resc)

# Multi-group comparison
resko1 <- richKEGG(gene1, kodata = hsako)
resko2 <- richKEGG(gene2, kodata = hsako)
res_cmp <- compareResult(list(S1 = resko1, S2 = resko2))
richCompareDot(res_cmp)

# Comparison heatmap
richHeatmap(list(GO = resgo, KEGG = resko), top = 50)
```

### UpSet plot

```r
gene_lists <- list(
  "Treatment A" = sample(unique(hsako$GeneID), 500),
  "Treatment B" = sample(unique(hsako$GeneID), 400),
  "Control"     = sample(unique(hsako$GeneID), 600)
)
richUpset(gene_lists)

# Custom colors and ordering
richUpset(gene_lists,
          mycol = c("dodgerblue", "goldenrod1", "seagreen3"),
          order.by = "degree", nintersects = 20)
```

## Working with results

```r
# Extract result table and detail table
result(resgo)
detail(resgo)

# Get genes for specific terms
getGenes(resgo, term = "GO:0006955")

# dplyr verbs work directly
library(dplyr)
resgo %>% filter(Padj < 0.01) %>% head()
resgo %>% select(Term, Pvalue, Padj)
resgo %>% arrange(Pvalue) %>% head(10)

# Batch enrichment across multiple gene lists
gene_lists <- list(GroupA = gene1, GroupB = gene2)
batch_res <- batchEnrich(gene_lists, annot = hsago)
```

## Function name aliases

All plotting functions use the `rich*` prefix as the primary name.
Old `gg*` names are kept as backward-compatible aliases:

| Primary name | Aliases |
|---|---|
| `richBar` | `ggbar` |
| `richDot` | `ggdot` |
| `richLollipop` | `gglollipop`, `ggLollipop` |
| `richGeneDot` | `gggenedot`, `ggGeneDot` |
| `richGeneHeat` | `gggeneheat`, `ggGeneHeat` |
| `richGeneBar` | `gggenebar`, `ggGeneBar` |
| `richCircle` | `ggcircbar`, `ggCircBar` |
| `richVolcano` | `ggvolcano`, `ggVolcano` |
| `richNES` | `ggNES`, `ggnes` |
| `richScatter` | `ggscatter`, `ggScatter` |
| `richECDF` | `ggecdf`, `ggECDF` |
| `richTermSim` | `ggtermsim`, `ggTermSim` |
| `richNetplot` | `ggnetplot` |
| `richNetwork` | `ggnetwork` |
| `richNetmap` | `ggnetmap` |
| `richUpset` | `ggupset` |
| `richHeatmap` | `ggheatmap` |
| `richClusterDot` | `ggcluster` |
| `richCompareDot` | `comparedot` |
| `richGSEAplot` | `ggGSEA` |
| `richGSEAcurve` | `plotGSEA` |

## Citation

If you use richR in your research, please cite:

```
Guo K and Hur J (2020). richR: Enrichment analysis for functional genomics.
R package version 0.1.2. https://github.com/guokai8/richR
DOI: 10.5281/zenodo.3675760
```

## Contact

For questions or bug reports please contact guokai8@gmail.com or open an
[issue on GitHub](https://github.com/guokai8/richR/issues).
