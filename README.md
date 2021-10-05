# PbPbRedux
TIMS Pb/Pb data reduction software

## Installation

To install `PbPbRedux` from GitHub, you need the `remotes` package. If
this is not installed on your system, then you can do so as follows:

```
install.packages('remotes')
```

Once `remotes` is installed, `PbPbRedux` can be installed as follows:

```
remotes::install_github('pvermees/PbPbRedux')
```

`PbPbRedux` depends on
[`IsoplotR`](http://github.com/pvermees/IsoplotR), so the latter
package is automatically loaded whenever `PbPbRedux` is.

## Example

Using the example input files in the [`inst`](inst) folder of this
GitHub repository:

```
library(PbPbRedux)

spikes <- read.csv('spikes.csv',header=TRUE)
settings('iratio','U238U235',137.786,0.015)

option <- 3

if (option==1){ # all samples use the same blank:
    samples <- read.csv('samples1.csv',header=TRUE)
    blanks <- read.csv('blanks1.csv',header=TRUE)
    tab <- process(samples,blanks,spikes)
} else if (option==2){ # each aliquot has its own blank:
    samples <- read.csv('samples2.csv',header=TRUE)
    blanks <- read.csv('blanks2.csv',header=TRUE)
    tab <- process(samples,blanks,spikes)
} else if (option==3){ # individual blanks with shared covariance matrix:
    samples <- read.csv('samples2.csv',header=TRUE)
    cblanks <- read.csv('blanks1.csv',header=TRUE)
    blanks <- read.csv('blanks2.csv',header=TRUE)
    tab <- process(samples,blanks,spikes,cblanks)
}

omit <- c(1:3,6,8) # aliquots to omit from the isochron

# plot the results using IsoplotR functions:
PbPb <- IsoplotR:::as.PbPb(tab[-omit,],format=2,ierr=4)
isochron(PbPb,exterr=FALSE)
```
