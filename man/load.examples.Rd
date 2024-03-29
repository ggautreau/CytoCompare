% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{load.examples}
\alias{load.examples}
\title{Retrieving of an example dataset of CytoCompare objects}
\usage{
load.examples(del.file = FALSE)
}
\arguments{
\item{del.file}{a logical specifying if an existing CytoCompareExample.rdata file can to be overwritten}
}
\value{
none
}
\description{
Downloads and loads an example dataset of CytoCompare objects constructed based on cytometry profiles obtained from healthy human bone marrow unstimulated or stimulated (PMID:21964415).

This example dataset consists on three cytometry profiles of healthy human bone marrow, unstimulated or stimulated by BCR or IL-7, measured using a mass cytometry panel of more than 30 cell markers. This panel has been designed to identify a large spectrum of immune cell types like monocytes, B, or CD4+ and CD8+ T cells. A SPADE analysis has been performed to identify cell clusters, that have been then manually labelled based on theirs profiles. SPADE cell clusters corresponding to 6 majors cell types have been extracted and a set of rectangle gates have been constructed based these cell types.

Once downloaded, the following objects will be available:\cr
* `bm_example.cells.b`, a `CELL` object containing the cell profiles of the B cell populations;\cr
* `bm_example.cells.mono`, a `CELL` object containing the cell profiles of the monocyte cell populations;\cr
* `bm_example.cells.tCD4naive`, a `CELL` object containing the cell profiles of the naive CD4+ T cell populations;\cr
* `bm_example.cells.tCD8naive`, a `CELL` object containing the cell profiles of the naive CD8+ T cell populations;\cr
* `bm_example.cells.tCD4mem`, a `CELL` object containing the cell profiles of the memory CD4+ T cell populations;\cr
* `bm_example.cells.tCD8mem`, a `CELL` object containing the cell profiles of the memory CD8+ T cell populations;\cr
* `bm_example.clusters`, a `CLUSTER` object containing the cell cluster profiles for all the different cell populations, identified by SPADE;\cr
* `bm_example.clusters.b`, a `CLUSTER` object containing the cell cluster profiles of the B cell populations, identified by SPADE;\cr
* `bm_example.clusters.mono`, a `CLUSTER` object containing the cell cluster profiles of the monocyte cell cluster profiles, identified by SPADE;\cr
* `bm_example.clusters.tCD4naive`, a `CLUSTER` object containing the cell cluster profiles of the naive CD4+ T cell populations, identified by SPADE;\cr
* `bm_example.clusters.tCD8naive`, a `CLUSTER` object containing the cell cluster profiles of the naive CD8+ T cell populations, identified by SPADE;\cr
* `bm_example.clusters.tCD4mem`, a `CLUSTER` object containing the cell cluster profiles of the memory CD4+ T cell populations, identified by SPADE;\cr
* `bm_example.clusters.tCD8mem`, a `CLUSTER` object containing the cell cluster profiles of the memory CD8+ T cell populations, identified by SPADE;\cr
* `bm_example.gates`, a `GATE` object containing the gate profiles constructed based on the six main cell populations identified by SPADE;\cr
* `bm_example.gates.b`, a `GATE` object containing the gate profiles constructed based on the B cell populations;\cr
* `bm_example.gates.mono`, a `GATE` object containing the gate profiles constructed based on the monocyte cell populations;\cr
* `bm_example.gates.tCD4naive`, a `GATE` object containing the gate profiles constructed based on the naive CD4T cell populations;\cr
* `bm_example.gates.tCD8naive`, a `GATE` object containing the gate profiles constructed based on the naive CD8T cell populations;\cr
* `bm_example.gates.tCD4mem`, a `GATE` object containing the gate profiles constructed based on the memory CD4T cell populations;\cr
* `bm_example.gates.tCD8mem`, a `GATE` object containing the gate profiles constructed based on the memory CD8T cell populations;\cr
* `bm_example.mweights`, a `MWEIGHTS` object containing cell markers that can be used in for comparison computations;\cr
* `bm_example.visne`, a list of three `CELL` objects containing the viSNE cell profiles for each biological sample.
}
\details{
This function downloads a CytoCompareExample.rdata file (from a public ftp server "ftp://ftp.cytocompare.org/public/rdata/") containing the different CytoCompare objects.
}

