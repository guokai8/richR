% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareResult.R
\name{comparedot}
\alias{comparedot}
\title{draw dotplot for multiple enrichment results}
\usage{
comparedot(
  x,
  pvalue = 0.05,
  low = "lightpink",
  high = "red",
  level = "Level3",
  pathway = NULL,
  alpha = 0.7,
  font.x = "bold",
  font.y = "bold",
  fontsize.x = 10,
  fontsize.y = 10,
  short = FALSE,
  padj = NULL,
  usePadj = TRUE,
  filename = NULL,
  width = 10,
  height = 8,
  include.all = FALSE
)
}
\arguments{
\item{x}{dataframe of enrichment result}

\item{pvalue}{cutoff value of pvalue (if padj set as NULL)}

\item{low}{low color}

\item{high}{high color}

\item{level}{the level for KEGG pathway}

\item{pathway}{the KEGG pathway name for the level}

\item{alpha}{transparency alpha}

\item{font.x}{font of x axis}

\item{font.y}{font of y axis}

\item{fontsize.x}{fontsize of x axis}

\item{fontsize.y}{fontsize of y axis}

\item{short}{automatic short name or not}

\item{padj}{cutoff value of p adjust value}

\item{usePadj}{use p adjust value as color or not (should use with padj)}

\item{filename}{figure output name}

\item{width}{figure width}

\item{height}{figure height}

\item{include.all}{include all richResults even empty}
}
\description{
draw dotplot for multiple enrichment results
}
\examples{
\dontrun{
hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
gene1 <- sample(unique(hsako$GeneID),1000)
gene2 <- sample(unique(hsako$GeneID),1000)
resko1 <-richKEGG(gene1,kodata = hsako)
resko2 <-richKEGG(gene2,kodata = hsako)
res<-compareResult(list(S1=resko1,S2=resko2))
comparedot(res,pvalue=0.05)
}
}
\author{
Kai Guo
}
