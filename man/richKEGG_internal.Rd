% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/richKEGG.R
\name{richKEGG_internal}
\alias{richKEGG_internal}
\title{KEGG Pathway Enrichment analysis function}
\usage{
richKEGG_internal(
  x,
  kodata,
  pvalue = 0.05,
  padj = NULL,
  ontology = "KEGG",
  organism = NULL,
  keytype = "SYMBOL",
  minSize = 2,
  maxSize = 500,
  minGSSize = 10,
  maxGSSize = 500,
  keepRich = TRUE,
  filename = NULL,
  padj.method = "BH",
  builtin = TRUE,
  sep = ","
)
}
\arguments{
\item{x}{vector contains gene names or dataframe with DEGs information}

\item{kodata}{GO annotation data}

\item{pvalue}{cutoff pvalue}

\item{padj}{cutoff p adjust value}

\item{ontology}{KEGG}

\item{organism}{organism}

\item{keytype}{keytype for input genes}

\item{minSize}{minimal number of genes included in significant terms}

\item{maxSize}{maximum number of genes included in significant terms}

\item{minGSSize}{minimal size of genes annotated by ontology term for testing.}

\item{maxGSSize}{maximal size of each geneSet for analyzing}

\item{keepRich}{keep terms with rich factor value equal 1 or not (default: TRUE)}

\item{filename}{output filename}

\item{padj.method}{pvalue adjust method(default:"BH")}

\item{sep}{character string used to separate the genes when concatenating}

\item{bulitin}{use KEGG bulit in KEGG annotation or not(set FALSE if you want use newest KEGG data)}
}
\description{
KEGG Pathway Enrichment analysis function
}
\author{
Kai Guo
}
