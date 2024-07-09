# Couch_SGE_BRCA2

This repository provides the VarCall model tool used in the analysis of high throughput CRISPR based saturated genome editing (SGE) data, one type of multiplexed assays of variant effects (MAVEs). It also includes related bioinformatics processing and statistical analysis code. The VarCall model mainly uses R and the R package rjags to perform the Bayesian two Gaussian components modeling.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Test and Run](#Test-and-run)
- [License](#license)

# Overview
``VarCall`` aims to provide a comprehensive statistical model for SGE data analysis, including the modeling of batch effects, and the prediction of variant effects. The package utilizes a Bayesian hierachical two Gaussian components modeling. The package should be able to run on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows) as long as R and related R Packages are installed.

# System Requirements
## Hardware requirements
`VarCall` package requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
### OS Requirements
This package is supported for *Linux*, but in principle should be able to run on all other platforms such as OS X and Windows. The package has been tested on the following systems:
+ Linux: Red Hat Enterprise Linux 8.8
+ Windows: Windows 10 Enterprise

### Package Dependencies
`VarCall` mainly depends on JAGS and the following R packages.

```
knitr, rjags, R2WinBUGS, R2jags, mgcv
```

# Installation Guide:
### Download or clone from github
```
git clone https://github.com/najiemayo/Couch_SGE_BRCA2
```

### Install packages
First install JAGS in the system following the link here: https://mcmc-jags.sourceforge.io/
Then within R, type the following:
```
install.packages(c("knitr", "rjags", "R2WinBUGS", "R2jags", "mgcv"))
```
This should be done within one miniute.

# Test and Run:
- To run a demo on a input file `combined.raw.20.tsv` using ~20% of the variants from the full BRCA2 exons data set:
  - `cd VarCall`
  - Start R by type R in the command line
  - Within R, type the following `library(knitr); knit("BRCA2mave24.ldaER.Rtex")`

- To run on the full data set:
  - Download the full data file `combined.raw.tsv` from GSE270424 and put it in the VarCall folder
  - Make sure the subfolder cache and figs are empty
  - Within file BRCA2mave24.ldaER.Rtex, change data set name by changing the line `% db<-"uvCounts"` to `% db<-"uvCountsFnl"`, and increase the MCMC iterations by changing 'mcmc.pars<-list(iter=10000, burn=5000, thin=10) to `mcmc.pars<-list(iter=150000, burn=50000, thin=10)`   
  - Start R by type R in the command line
  - Within R, type the following `library(knitr); knit("BRCA2mave24.ldaER.Rtex")`
    
- Specified prior:
  - Currently the prior is set to a mean value of 0.2 using a beta distribution Beta(2, 8). The change the prior, modifiy the line `beta.a<-2.0` and `beta.b<-8.0`. 

- Expected output and running time.
The output files will be `BRCA2mave24.ldaER.tex`, `MAVEpostProbs.csv` and several pdf plots. File `BRCA2mave24.ldaER.tex` and the plots can be further compiled into a pdf file. File `MAVEpostProbs.csv` is the main output file for predicted probabilities of being pathogenic based on the training labels from file `variant_type_for_train.csv`. The following shows the main columns of output:
  - PrDel: the probability of being pathogenic
  - lPostOdds: the log posterior odds
  - logBF: the log Bayes factor
  - eta: the estimated effect size
  - eta.ll, eta.ul: 95% lower bound and upper bound of eta

# License

This project is covered under the **MIT License**.


