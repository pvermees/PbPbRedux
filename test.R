setwd('~/Documents/Programming/R/PbPb_Redux/')
source('PbPb_Redux.R')
library(IsoplotR)

samples <- read.csv('samples.csv',header=TRUE)
blanks <- read.csv('blanks2.csv',header=TRUE)
spikes <- read.csv('spikes.csv',header=TRUE)

if (FALSE){
    i <- 3
    ablk <- avgblank(blanks,spikes,spk=samples[i,'spk']) # with covariance matrix
    lsmp <- getsample(i=i,samples,spikes) # without covariance matrix
    E <- getE(i,samples,ablk,spikes)
    pinit <- init(i=i,lsmp=lsmp,lblk=ablk$lblk)
    fit <- optim(pinit,fn=LL,i=i,lsmp=lsmp,lblk=ablk$lblk,E=E,hessian=TRUE)
    covmat <- solve(fit$hessian)
}
if (TRUE){
    tab <- process(samples[-2,],blanks,spikes)
    PbPb <- IsoplotR:::as.PbPb(tab,format=2)
    settings('iratio','U238U235',137.786,0.015)
    isochron(PbPb,exterr=FALSE)
}
