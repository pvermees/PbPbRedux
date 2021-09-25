# PbPbRedux
TIMS Pb/Pb data reduction software

Example:

```
# load the PbPbRedux code:
source('PbPbRedux.R')

# load the input files:
samples <- read.csv('samples.csv',header=TRUE)
blanks <- read.csv('blanks2.csv',header=TRUE)
spikes <- read.csv('spikes.csv',header=TRUE)

# process the data, omitting the 2nd aliquot:
tab <- process(samples[-2,],blanks[,-2],spikes)

# plot the isochron
PbPb <- IsoplotR:::as.PbPb(tab,format=2)
IsoplotR::isochron(PbPb)
```
