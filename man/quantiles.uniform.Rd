% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{quantiles.uniform}
\alias{quantiles.uniform}
\title{Quantiles of a uniform distribution}
\usage{
quantiles.uniform(bounds, value)
}
\arguments{
\item{bounds}{a numeric vector indicating the support (lower and upper bounds) of the uniform distribution}

\item{value}{a numeric value specifying the quantile to compute}
}
\value{
a numeric value of the quantiles value
}
\description{
Provides the quantiles of a uniform distribution
}
\details{
This function is used internally when comparing the marker expression ranges of gate profiles with the default comparison approach. This function can also be used by users willing to define their own statistical functions for comparing cell, cell cluster or gate profiles.
}

