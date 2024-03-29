% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{quantiles.density}
\alias{quantiles.density}
\title{Quantiles of a DENSITY object}
\usage{
quantiles.density(density, values)
}
\arguments{
\item{density}{a DENSITY object}

\item{values}{a numeric value or a numeric vector specifying the quantiles to compute}
}
\value{
a numeric value of the quantiles values
}
\description{
Provides the quantiles of a DENSITY object, at specific values.
}
\details{
This function is used internally when comparing the expression densities of cell markers with the proposed comparison approach. This function can also be used by users willing to define their own statistical functions for comparing cell, cell cluster or gate profiles.
}

