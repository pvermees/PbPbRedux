setwd('~/Documents/Programming/R/PbPbRedux/')
source('PbPbRedux.R')
library(IsoplotR)

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

# aliquots to omit from the isochron
omit <- c(1:3,6,8)
PbPb <- IsoplotR:::as.PbPb(tab[-omit,],format=2)
isochron(PbPb,exterr=FALSE)
