% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{ReadPolyApipe}
\alias{ReadPolyApipe}
\title{The function reads in output from polyApipe and returns a count matrix}
\usage{
ReadPolyApipe(
  counts.file,
  peaks.file = NULL,
  sep = c(":", ",", ":"),
  filter.chromosomes = TRUE,
  min.features = NULL,
  min.cells = NULL,
  verbose = TRUE
)
}
\arguments{
\item{counts.file}{file containing the counts from polyApipe in .tsv format}

\item{peaks.file}{gff file that requires the seqid, start, end, and strand of each peak}

\item{sep}{Separators}

\item{filter.chromosomes}{If TRUE (default), only include main chromosomes.}

\item{min.features}{Cell must have at least this many feature to be included.}

\item{min.cells}{A feature must be present in at least this manhy cells.}

\item{verbose}{Print output.}
}
\value{
Returns count matrix with peaks in the form "chromosome:start,end:strand"
such as "10:100000560,100000859,100000560,100000859:+"
}
\description{
The function reads in output from polyApipe and returns a count matrix
}
\concept{preprocessing}
