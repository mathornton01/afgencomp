% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afgencomp.R
\name{countMersCertainFirstFive}
\alias{countMersCertainFirstFive}
\title{Compute Vector of counts for first five kinds of kmers}
\usage{
countMersCertainFirstFive(genomicString, freq = FALSE)
}
\arguments{
\item{genomicString}{The character string version of the genomic signal which
the user wishes to produce the kmer vector for.}

\item{freq}{Whether the user would like a vector of raw counts back or if
a vector of relative frequencies.}
}
\value{
a 1364 length vector containing counts/frequencies for first five mers
}
\description{
Generates a vector of counts (or frequencies) for the first five k-mers (ie.
1-mers, 2-mers, ..., 5-mers) and returns the vector concatenated together.
}
\examples{
countMersCertainFirstFive(sarscvmay[1,]$sequences,FALSE)
}
