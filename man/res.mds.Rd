% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/res_display.R
\name{res.mds}
\alias{res.mds}
\title{Create a Multidimensional scaling (MDS) representation of a RES object}
\usage{
res.mds(res, filename = "res.html", cols = NULL, sizes = NULL,
  svgsize = 1000)
}
\arguments{
\item{res}{a RES object}

\item{filename}{a file location where to save the objects}

\item{cols}{a character vector specifying the colours of the node in the SVG representation}

\item{sizes}{a numeric vector specifying the sizes of the nodes in pixels in the SVG representation}

\item{svgsize}{a numeric value specifying the size of the SVG representation in pixels}
}
\value{
none
}
\description{
Creates a Multidimensional scaling (MDS) representation of the comparison results. In such MDS representation each cytometry profile is represented by a dot in a two-dimensional space and the distances between the nodes are proportional to the similarity measures between the profiles. The Kruskal Stress displayed at the left bottom of the representation quantifies the quality of the representation as the percentage of information lost in the dimensionality reduction process.
}
\details{
The representation is provided as a HTML file via a Scalable Vector Graphics (SVG) element created with the D3.js library.
}

