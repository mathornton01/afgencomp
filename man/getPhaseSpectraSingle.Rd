% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/afgencomp.R
\name{getPhaseSpectraSingle}
\alias{getPhaseSpectraSingle}
\title{Get the Fourier Phase Spectrum for a single encoded Genomic Signal}
\usage{
getPhaseSpectraSingle(encodedSignal)
}
\arguments{
\item{encodedSignal}{is a genomic signal of interest, for which the average
power spectrum will be computed and returned to the user, you may encode
your genomic character strings by using the encodeGenomes or encodeGenome
function in this package.}
}
\value{
A vector of values indicating the average phase spectral density
(according to the Fourier Transform) for the encoded genomes four, or two
constituent signals.
}
\description{
This function produces the phase spectrum of the Fourier transform for a
single genomic signal that has been encoded using either the 2D or 4D
representation, the function will produce an error if it is not supplied with
a matrix of values that has a number of rows equal to 2 or 4. (AM)
}
\examples{
MTHFR100 <- "TGGCCAGGTATAGTGGCTCATACCTGTAATCCCAGCACTCTGGGAGACCGAAGCAGTATCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATG";
encMTHFR100 <- encodeGenome(MTHFR100);
psMTHFR100 <- getPhaseSpectraSingle(encMTHFR100);
plot(psMTHFR100, type='l',xlab='Frequency/Sequency',
     ylab='Phase Spectral Density', main="Phase Spectrum of First 100 nucleotides of MTHFR");
}
