% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{utilities-methods}
\alias{utilities-methods}
\alias{head,Annot-method}
\alias{head,richResult-method}
\alias{head,GSEAResult-method}
\alias{tail,Annot-method}
\alias{tail,richResult-method}
\alias{tail,GSEAResult-method}
\alias{dim,Annot-method}
\alias{dim,richResult-method}
\alias{dim,GSEAResult-method}
\title{S4 methods for head, tail, dim on Annot, richResult, GSEAResult}
\usage{
\S4method{head}{Annot}(x, n = 6L, ...)

\S4method{head}{richResult}(x, n = 6L, ...)

\S4method{head}{GSEAResult}(x, n = 6L, ...)

\S4method{tail}{Annot}(x, n = 6L, ...)

\S4method{tail}{richResult}(x, n = 6L, ...)

\S4method{tail}{GSEAResult}(x, n = 6L, ...)

\S4method{dim}{Annot}(x)

\S4method{dim}{richResult}(x)

\S4method{dim}{GSEAResult}(x)
}
\arguments{
\item{x}{An object of class \code{Annot}, \code{richResult}, or \code{GSEAResult}}

\item{n}{Number of rows to display (for \code{head} or \code{tail})}

\item{...}{Further arguments passed to the underlying data.frame method}
}
\description{
Defines how to call \code{head}, \code{tail}, and \code{dim} for these S4 classes.
}
