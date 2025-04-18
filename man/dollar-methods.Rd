% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{dollar-methods}
\alias{dollar-methods}
\alias{$,Annot-method}
\alias{$-Annot}
\alias{$,richResult-method}
\alias{$-richResult}
\alias{$,GSEAResult-method}
\alias{$-GSEAResult}
\title{S4 `$` accessors for Annot and richResult}
\usage{
\S4method{$}{Annot}(x, name)

\S4method{$}{richResult}(x, name)

\S4method{$}{GSEAResult}(x, name)
}
\arguments{
\item{x}{An S4 object (Annot or richResult)}

\item{name}{column name to extract}
}
\description{
Provide `obj$name` access that extracts from the relevant slot columns.
}
