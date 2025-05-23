% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/richKEGG.R
\name{richKEGG,ANY,data.frame-method}
\alias{richKEGG,ANY,data.frame-method}
\title{KEGG Pathway Enrichment analysis function}
\usage{
\S4method{richKEGG}{ANY,data.frame}(
  x,
  kodata,
  pvalue = 0.05,
  padj = NULL,
  organism = NULL,
  ontology = "KEGG",
  keytype = NULL,
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

\item{pvalue}{cutoff pvalue}

\item{padj}{cutoff p adjust value}

\item{organism}{organism}

\item{ontology}{KEGG}

\item{keytype}{keytype for input genes}

\item{minSize}{minimal number of genes included in significant terms}

\item{maxSize}{maximum number of genes included in significant terms}

\item{minGSSize}{minimal size of genes annotated by ontology term for testing.}

\item{maxGSSize}{maximal size of each geneset for analyzing}

\item{keepRich}{keep terms with rich factor value equal 1 or not (default: TRUE)}

\item{filename}{output filename}

\item{padj.method}{pvalue adjust method(default:"BH")}

\item{sep}{character string used to separate the genes when concatenating}

\item{bulitin}{use KEGG bulit in KEGG annotation or not(set FALSE if you want use newest KEGG data)}
}
\description{
KEGG Pathway Enrichment analysis function
}
\examples{
\dontrun{
  hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
  hsako<-as.data.frame(hsako)
  gene=sample(unique(hsako$GeneID),1000)
  res<-richKEGG(gene,kodata = hsako)
}
}
\author{
Kai Guo
}
