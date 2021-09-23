# PbPb_Redux
TIMS Pb/Pb data reduction software

Example:

```
# load the PbPb_Redux code:
source('PbPb_Redux.R')

# load the input files:
samples <- read.csv('samples.csv',header=TRUE)
blanks <- read.csv('blanksGujba.csv',header=TRUE)
spikes <- read.csv('spikes.csv',header=TRUE)

# process the data, omitting the 2nd aliquot:
tab <- process(samples[-2,],blanks,spikes)

# plot the isochron
PbPb <- IsoplotR:::as.PbPb(tab,format=2)
IsoplotR::isochron(PbPb)
```
