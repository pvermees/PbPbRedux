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

## Examples

Using the example input files in the [`inst`](inst) folder of this
GitHub repository:

### Example 1: all samples use the same blank

```
# load the package:
library(PbPbRedux)

# change the default 238U/235U ratio to a meteoritic value:
settings('iratio','U238U235',137.786,0.015)

# load the input files:
spikes <- read.csv('spikes.csv',header=TRUE)
samples <- read.csv('samples1.csv',header=TRUE)
blanks <- read.csv('blanks1.csv',header=TRUE)

# process the data:
tab <- process(samples,blanks,spikes,ierr=4)

# plot the results using IsoplotR functions:
PbPb <- IsoplotR:::as.PbPb(tab,format=2,ierr=4)
isochron(PbPb,exterr=FALSE)
```

### Example 2: each aliquot has its own blank

This is the same code as the previous example, but using different
input files. See the [`inst`](inst) folder or details:

```
library(PbPbRedux)
settings('iratio','U238U235',137.786,0.015)
spikes <- read.csv('spikes.csv',header=TRUE)
samples <- read.csv('samples2.csv',header=TRUE)
blanks <- read.csv('blanks2.csv',header=TRUE)
tab <- process(samples,blanks,spikes)
PbPb <- IsoplotR:::as.PbPb(tab,format=2,ierr=4)
isochron(PbPb,exterr=FALSE)
```

### Example 3: individual blanks with shared covariance matrix

Now we load an extra file of blanks, whose covariance matrix is used
as a measure of dispersion for the individual blank measurements:

```
library(PbPbRedux)
settings('iratio','U238U235',137.786,0.015)
spikes <- read.csv('spikes.csv',header=TRUE)
samples <- read.csv('samples2.csv',header=TRUE)
blanks <- read.csv('blanks2.csv',header=TRUE)

# load 15 replicate blank values:
cblanks <- read.csv('blanks1.csv',header=TRUE)

# process the data and plot the results:
tab <- process(samples,blanks,spikes,cblanks)
PbPb <- IsoplotR:::as.PbPb(tab,format=2,ierr=4)
isochron(PbPb,exterr=FALSE)
```

### Omitting aliquots from the isochron

Individual aliquots can be removed from the isochron using standard
`R` syntax. Following up from any of the above examples:

```
omit <- c(1:3,6,8) # aliquots to omit from the isochron
PbPb <- IsoplotR:::as.PbPb(tab[-omit,],format=2,ierr=4)
isochron(PbPb,exterr=FALSE)
```

