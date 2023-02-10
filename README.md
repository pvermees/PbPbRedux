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

## Examples

Using the example input files in the [`inst`](inst) folder of this
GitHub repository:

### Example 1: all samples use the same blank

```
# load the package:
library(PbPbRedux)

# load the input files:
smp <- read.csv('samples1.csv',header=TRUE)
blk <- read.csv('blanks1.csv',header=TRUE)
spk <- read.csv('spikes.csv',header=TRUE)

# process the data:
tab <- process(samples=smp,blanks=blk,spikes=spk)
```

### Example 2: each aliquot has its own blank

This is the same code as in the previous example, but using different
input files. See the [`inst`](inst) folder for details:

```
library(PbPbRedux)
smp <- read.csv('samples2.csv',header=TRUE)
blk <- read.csv('blanks2.csv',header=TRUE)
spk <- read.csv('spikes.csv',header=TRUE)
tab <- process(samples=smp,blanks=blk,spikes=spk)
```

### Example 3: individual blanks with shared covariance matrix

Now we load an extra file of blanks, whose covariance matrix is used
as a measure of dispersion for the individual blank measurements:

```
library(PbPbRedux)
smp <- read.csv('samples2.csv',header=TRUE)
blk <- read.csv('blanks2.csv',header=TRUE)
spk <- read.csv('spikes.csv',header=TRUE)

# load 15 replicate blank values:
cblk <- read.csv('blanks1.csv',header=TRUE)

tab <- process(samples=smp,blanks=blk,spikes=spk,cblanks=cblk)
```

To view the documentation of the `process` function:

```
?process
```

### Saving and plotting the results:

The output table can be saved to a comma separated variable
file using the basic `write.csv` function:

```
write.csv(tab,file='results.csv',row.names=FALSE)
```

You can either open this file in a spreadsheet application for further
processing, or you can continue working in `R`.  For example, plotting
the results as an isochron using the
[`IsoplotR`](http://github.com/pvermees/IsoplotR) package:

```
library(IsoplotR)

# change the default 238U/235U ratio to a meteoritic value:
settings('iratio','U238U235',137.786,0.015)

# using IsoplotR's undocumented 'as.PbPb' function to convert
# the output table into a PbPb data object, where ierr=4 specifies
# that errors are provided as 2se % values:
PbPb <- IsoplotR:::as.PbPb(tab,format=2,ierr=4)

# plot as an inverse isochron (type ?isochron for additional options):
isochron(PbPb,exterr=FALSE)
```

Individual aliquots can be removed from the isochron using standard
`R` syntax:

```
omit <- c(1:3,6,8) # aliquots to omit from the isochron
PbPb <- IsoplotR:::as.PbPb(tab[-omit,],format=2,ierr=4)
isochron(PbPb,exterr=FALSE)
```