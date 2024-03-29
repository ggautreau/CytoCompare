% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{methods}
\name{extract}
\alias{[,CELL,ANY,ANY-method}
\alias{[,CLUSTER,ANY,ANY-method}
\alias{[,GATE,ANY,ANY-method}
\alias{[,MWEIGHTS,ANY,ANY-method}
\alias{[,RES,ANY,ANY-method}
\alias{extract}
\title{Extraction of subsets of data from CytoCompare objects}
\usage{
\S4method{[}{CELL,ANY,ANY}(x, i, j)

\S4method{[}{CLUSTER,ANY,ANY}(x, i, j)

\S4method{[}{GATE,ANY,ANY}(x, i, j)

\S4method{[}{MWEIGHTS,ANY,ANY}(x, i)

\S4method{[}{RES,ANY,ANY}(x, i)
}
\arguments{
\item{x}{a CELL, CLUSTER, GATE, MWEIGHTS or RES object}

\item{i}{a numeric, logical or character vector}

\item{j}{a numeric, logical or character vector}
}
\value{
a S4 object of class CELL, CLUSTER, GATE, MWEIGHTS or RES
}
\description{
Extracts subsets of CELL, CLUSTER, GATE, MWEIGHTS or RES object.
}
\details{
For cytometry objects (CELL, CLUSTER, or GATE objects), the parameter i represents a vector of profiles to extract and the parameter j represents a vector of markers to extract.

For MWEIGHTS objects, the parameter i represents a vector of markers to extract.

For RES objects, the parameter i represents a vector of comparisons to extract.
}

