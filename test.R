setwd('~/Documents/Programming/R/PbPb_Redux/')
source('PbPb_Redux.R')

samples <- read.csv('samples.csv',header=TRUE)
blanks <- read.csv('blanksGujba.csv',header=TRUE)
spikes <- read.csv('spikes.csv',header=TRUE)

if (FALSE){
    i <- 3
    ablk <- avgblank(blanks,spikes,spk=samples[i,'spk']) # with covariance matrix
    lsmp <- getsample(i=i,samples,spikes) # without covariance matrix
    E <- getE(i,samples,ablk,spikes)
    pinit <- init(i=i,lsmp=lsmp,lblk=ablk$lblk)
    fit <- optim(pinit,fn=LL,i=i,lsmp=lsmp,lblk=ablk$lblk,E=E,hessian=TRUE)
    covmat <- solve(fit$hessian)
} else {
    tab <- process(samples[-c(2,14),],blanks,spikes)
    PbPb <- IsoplotR:::as.PbPb(tab,format=2)
    IsoplotR::isochron(PbPb)
}
