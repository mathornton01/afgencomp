# `Afgencomp`: An R-Package for Alignment-Free Genetic Comparison (Previously `YinGenomicDFTDistance`). 

This repository contains a software distribution of an R-Package, Afgencomp, which contains utility functions for comparing string representations of 
genomic sequences, without having to perform multiple sequence alignment (MSA, an often timely and expensive step).  

## Package Contents
The package contains documented functions which allow the user to: 

* Determine the optimal (balance preserving) 2-D encoding strategy for a genomic sequence (`getOptimal2DStrategy()`).
* Encode (a) genomic sequence(s) using a 4-D Binary representation, or 2-D 3-level representation. (`encodeGenome()`, `encodeGenomes()`).
* Calculate the Fourier power spectra for (a) sequence(s) (`getPowerSpectraSingle()`,`getPowerSpectraEnsemble()`). 
* *Added Post-Dissertation Defense:* Calculate the Fourier phase spectra for (a) sequence(s) (`getPhaseSpectraSingle()`,`getPhaseSpectraEnsemble()`). 
* Evenly scale (a) series (of length $n$ < $m$) to a specified length $m$ (`evenlyScaleSingle()`,`evenlyScaleEnsemble()`). 
* Directly compute the euclidean distance in power spectral density matrix for all genomes in list (utility function - `getPowerSpectraDistances()`).
* Determine k-mer counts for small values of $k$ (`countKmersCertain()`).
* Compute Vector of counts for first five kinds of kmers (utility functions - `countMersCertainFirstFive()`, `countFFMersEnsemble()`).
* Filter components of Power/Phase spectra with low variance across genomes in sample (`getMinimalVarianceFilter()`).
* One test function for checking the package integrity (`tstForMultiGenomes()`). 


