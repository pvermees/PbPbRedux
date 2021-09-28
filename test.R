setwd('~/Documents/Programming/R/PbPbRedux/')
source('PbPbRedux.R')
library(IsoplotR)

samples <- read.csv('samples.csv',header=TRUE)
blanks <- read.csv('blanks2.csv',header=TRUE)
spikes <- read.csv('spikes.csv',header=TRUE)

tab <- process(samples[-2,],blanks[-2,],spikes)
PbPb <- IsoplotR:::as.PbPb(tab,format=2)
settings('iratio','U238U235',137.786,0.015)
isochron(PbPb,exterr=FALSE,show.numbers=TRUE)
