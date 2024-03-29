% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imports.R
\name{import.SPADE}
\alias{import.SPADE}
\title{Importation of cell cluster profiles from SPADE results}
\usage{
import.SPADE(path, dictionary = NULL, exclude = NULL, trans = "arcsinh",
  trans.para = list(arcsinh.scale = 5), bin.width = 0.05,
  extract.folder = NULL, extract.folder.del = FALSE, zip = FALSE)
}
\arguments{
\item{path}{a character indicating the location to a zip or a folder containing the SPADE results}

\item{dictionary}{a two-column data.frame providing the correspondence between the original marker names (first column) and the new marker names (second column)}

\item{exclude}{a character vector containing the marker names to be excluded in the import procedure}

\item{trans}{a character specifying the name of a transformation function to apply on the marker expression intensities. Possible functions are "arcsinh" for arc sin hyperbolic transformation (default), "log" for logarithmic transformation, or "none" for no transformation}

\item{trans.para}{a character named list containing parameters for the transformation. Please refer to the details section for more details}

\item{bin.width}{a numeric value indicating the width of the bins for the marker expression densities computations}

\item{extract.folder}{a folder path for extracting the SPADE zip archive (temporary folder by default)}

\item{extract.folder.del}{a logical value indicating if the extracted SPADE results should be removed after the extraction}

\item{zip}{a logical value that specify if the path specifies a zip file}
}
\value{
a S4 object of class CLUSTER
}
\description{
Imports one or several cell cluster profiles identified by the SPADE algorithm into a CLUSTER object.
}
\details{
SPADE is a popular visualization and analysis algorithm that identifies clusters of cells having similar expression profiles for selected markers using an agglomerative hierarchical clustering-based algorithm combined with a density-based down-sampling procedure (PMID:21964415). Given a set of FCS files (usually one file per sample), SPADE identifies cell clusters based on the whole dataset and provides then for each sample the amount of cells present within each cluster.
}

