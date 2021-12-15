# `Afgencomp`: An R-Package for Alignment-Free Genetic Comparison (Previously `YinGenomicDFTDistances`). 

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

If you are looking to install
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

```bash
git clone https://github.com/mathornton01/afgencomp.git 
```

Once the package has been downloaded it can be installed for R by navigating to 
the `afgencomp` directory from within R and using the `install.packages()`
function. 

```R
afgencomp.pkg.dir <- "full/path/to/directory/goes/here";
install.packages(afgencomp.pkg.dir, repos=NULL, type="source");
````

Or if you would instead like to choose the file using the file-manager, you 
can do so by running: 

```R
install.packages(file.choose(), repos=NULL)
```

then selecting the top-level directory for the package source.  That is the 
directory which contains the `R/`, `data/`, `man/`, and `vignettes/` folders. 

## Install The Software from Its Official Github Repository Directly

The `devtools` library in R allows for developers to quickly and easily share 
there packages with R-users via Github.  the `install_github` function of the 
`devtools` package.  Be sure to specify that you would like for the package 
vignette (this document) to be constructed when you run this, so that the 
vignette is available via the `browseVignettes()` or `??` functions. 

```R
library(devtools);
```

```R
devtools::install_github("https://github.com/mathornton01/afgencomp.git",build_vignettes = TRUE);
```

# References 

Much of this package was intended to implement the genomic comparison strategy implemented and used by Yin and Yau since 
they published it in 2014.  A selection of the papers describing the techniques for power spectral distance calculation and
scaling that are included in the package are provided: 

* An improved model for whole genome phylogenetic analysis by Fourier transform, *Yin, Changchuan and Yau, Stephen S-T*, **Journal of theoretical biology vol. 382, pp. 99-110,** 2015. Elsevier.
* A measure of DNA sequence similarity by Fourier Transform with applications on hierarchical clustering, *Yin, Changchuan and Chen, Ying and Yau, Stephen S-T*, **Journal of theoretical biology vol. 359, pp.18-28,** 2014. Elsevier.
* A new method to cluster DNA sequences using Fourier power spectrum, *Hoang, Tung and Yin, Changchuan and Zheng, Hui and Yu, Chenglong and He, Rong Lucy and Yau, Stephen S-T*, **Journal of theoretical biology vol. 372, pp. 135-145,** 2015. Elsevier. 
*  CPhylogenetic analysis of DNA sequences or genomes by Fourier transform, *Changchuan Yin* (2021). **MATLAB Central File Exchange. Retrieved December 15, 2021.** [link](https://www.mathworks.com/matlabcentral/fileexchange/52072-phylogenetic-analysis-of-dna-sequences-or-genomes-by-fourier-transform),  


