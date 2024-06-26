Library dependencies: knitr, rjags, R2WinBUGS, R2jags, mgcv

The script itself is the file `BRCA2mave24.ldaER.Rtex'.  To compile this, get into the folder VarCallMAVEforGitHub, start R, then run the following in R:

library(knitr)
knit("BRCA2mave24.ldaER.Rtex")    

Run this on a machine with a fair amount of RAM memory and a reasonably fast CPU clock.  It will take 6 to 12 hours, I believe. You can check the code by changing some of the user-specified parameters at the beginning of the script:  if you set `subst<-TRUE' it will run using a random subset of 20% of the data and if you choose the first `mcmc.pars' choice (with `iters=10000') you can check that it works on your system likely less than 1 hour. You can then clear the cache, reset these parameters then run the model using the same settings used for the paper.