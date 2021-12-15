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

# Package Installation

This vignette is distributed and maintained along with the software in `afgencomp`. 
However, if you happen to find your way to this vignette, and are looking to install 
the `afgencomp` R-package you may do so via one of several routes: 

* Install the Software from A Source Distribution using `install.packages()`
* Install the Software from the Github Repository using `devtools::install_github()`
* (Coming Soon) Install the Software from CRAN using `install.packages()`

## Installing The Software from Source 

The source distribution of the software is written in R, and hosted on Github. 
The source code can be downloaded using the git software to clone the `afgencomp`
repository.  This can be accomplished by simply typing the following `git` 
command in a terminal window where you would like to keep the package files 
or by downloading the package manually using a web-browser and navigating to the 
[github page](https://github.com/mathornton01/afgencomp). 

`
bash> git clone https://github.com/mathornton01/afgencomp.git 
`

Once the package has been downloaded it can be installed for R by navigating to 
the `afgencomp` directory from within R and using the `install.packages()`
function. 

`
r> afgencomp.pkg.dir <- "full/path/to/directory/goes/here";
r> install.packages(afgencomp.pkg.dir, repos=NULL, type="source");
`

Or if you would instead like to choose the file using the file-manager, you 
can do so by running: 

`
r> install.packages(file.choose(), repos=NULL)
`

then selecting the top-level directory for the package source.  That is the 
directory which contains the `R/`, `data/`, `man/`, and `vignettes/` folders. 

## Install The Software from Its Official Github Repository Directly

The `devtools` library in R allows for developers to quickly and easily share 
there packages with R-users via Github.  the `install_github` function of the 
`devtools` package.  Be sure to specify that you would like for the package 
vignette (this document) to be constructed when you run this, so that the 
vignette is available via the `browseVignettes()` or `??` functions. 
`
r> library(devtools);
r> devtools::install_github("https://github.com/mathornton01/afgencomp.git",build_vignettes = TRUE);
`

