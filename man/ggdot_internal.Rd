% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggdot.R
\name{ggdot_internal}
\alias{ggdot_internal}
\title{Dotplot for enrichment results}
\usage{
ggdot_internal(
  object,
  top = 50,
  pvalue = 0.05,
  order = FALSE,
  low = "lightpink",
  high = "red",
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
  height = 8
)
}
\arguments{
\item{object}{richResult object}

\item{top}{number of terms you want to display,}

\item{pvalue}{cutoff value of pvalue (if padj set as NULL)}

\item{low}{low color}

\item{high}{high color}

\item{alpha}{transparency alpha}

\item{font.x}{font of x axis}

\item{fontsize.x}{fontsize of x axis}

\item{fontsize.y}{fontsize of y axis}

\item{short}{automatic short name or not}

\item{padj}{cutoff value of p adjust value}

\item{usePadj}{use p adjust value as color or not (should use with padj)}

\item{filename}{figure output name}

\item{width}{figure width}

\item{height}{figure height}

\item{font}{y font of y axis}

\item{font.size}{font size for xlim or ylim}
}
\description{
Dotplot for enrichment results
}
