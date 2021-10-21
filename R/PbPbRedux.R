# convert input errors to relative uncertainties
ierradj <- function(dat,rcol,ierr=1){
    errcol <- rcol + 1
    out <- dat
    if (ierr==1){
        out[,errcol] <- dat[,errcol]/dat[,rcol]
    } else if (ierr==2){
        out[,errcol] <- 0.5*dat[,errcol]/dat[,rcol]
    } else if (ierr==3){
        out[,errcol] <- dat[,errcol]/100
    } else if (ierr==4){
        out[,errcol] <- dat[,errcol]/200
    } else {
        stop('Invalid ierr value.')
    }
    out
}

# convert standard errors to specified output format
oerradj <- function(dat,rcol,oerr=1){
    errcol <- rcol + 1
    out <- dat
    if (oerr==1){
        # do nothing
    } else if (oerr==2){
        out[,errcol] <- 2*dat[,errcol]
    } else if (oerr==3){
        out[,errcol] <- 100*dat[,errcol]/dat[,rcol]
    } else if (oerr==4){
        out[,errcol] <- 200*dat[,errcol]/dat[,rcol]
    } else {
        stop('Invalid oerr value.')
    }
    out
}

beta <- function(mx){
    (log(205)-log(mx))/(log(205)-log(202))
}

lxstar <- function(lM,lr45,lr65,lr75,lr85,lr25,lr25t){
    if (missing(lr25t)) fract <- -lr25
    else fract <- lr25t-lr25
    l4 <- lr45 + lM + beta(204)*fract
    l6 <- lr65 + lM + beta(206)*fract
    l7 <- lr75 + lM + beta(207)*fract
    l8 <- lr85 + lM + beta(208)*fract
    out <- cbind(l4,l6,l7,l8)
    colnames(out) <- c('l4','l6','l7','l8')
    out
}

avgblank <- function(i,blanks,blk,spikes,spk,conc=TRUE){
    if (missing(blk)) blk <- blanks[i,'name']
    if (missing(spk)) spk <- blanks[i,'spk']
    if (missing(i)) i <- which(blanks[,'spk']%in%spk & blanks[,'name']%in%blk)
    if (length(i)<1) stop('Missing blank data.')
    lbdat <- log(blanks[i,c('mgspk','r52','r54','r56','r57','r58'),drop=FALSE])
    ispk <- which(spk%in%spikes[,'name'])
    lblk <- lxstar(lM=lbdat[,'mgspk'] +
                       ifelse(conc,log(spikes[ispk,'pmg205']),0),
                   lr45=-lbdat[,'r54'], # flip sign
                   lr65=-lbdat[,'r56'], # flip sign
                   lr75=-lbdat[,'r57'], # flip sign
                   lr85=-lbdat[,'r58'], # flip sign
                   lr25=-lbdat[,'r52'], # flip sign
                   lr25t=-log(spikes[spk,'r52'])) # flip sign
    lblkprime <- lxstar(lM=lbdat[,'mgspk'] +
                            ifelse(conc,log(spikes[ispk,'pmg205']),0),
                        lr45=-lbdat[,'r54'], # flip sign
                        lr65=-lbdat[,'r56'], # flip sign
                        lr75=-lbdat[,'r57'], # flip sign
                        lr85=-lbdat[,'r58'], # flip sign
                        lr25=-lbdat[,'r52']) # flip sign
    mlblk <- colMeans(lblkprime)
    if (length(i)<2){
        E <- diag(blanks[i,c('mgspkerr','err52','err54',
                             'err56','err57','err58')])^2
        J <- matrix(0,4,6)
        rownames(J) <- c('l4','l6','l7','l8')
        J[1:4,1] <- 1
        J[1,2] <- -beta(204)
        J[2,2] <- -beta(206)
        J[3,2] <- -beta(207)
        J[4,2] <- -beta(208)
        J[1,3] <- 1
        J[2,4] <- 1
        J[3,5] <- 1
        J[4,6] <- 1
        covlblk <- J %*% E %*% t(J)
    } else {
        covlblk <- cov(lblkprime)
    }
    lspk <- -log(spikes[ispk,'r52'])[1]
    errlspk <- 0 # placeholder for spikes[spk[1],'err52']/200
    out <- list(lblk=mlblk,covlblk=covlblk,lspk=lspk,errlspk=errlspk)
    if (conc){
        pgPb <- sweep(lblk,MARGIN=2,FUN='+',log(c(204,206,207,208)))
        pgPbt <- rowSums(exp(pgPb))
        out$Pb <- exp(mean(log(pgPbt)))
        if (length(i)<2){
            J <- exp(pgPb)
            out$relerrPb <- sqrt(J %*% covlblk %*% t(J))/pgPbt
        } else {
            out$relerrPb <- sd(log(pgPbt))
        }
    }
    out
}

# get x/5 data from a particular aliquot:
getaliquot <- function(i,samples,spikes,conc=TRUE){
    inames <- c('mgspk','r74','r64','r76','r65','r86','r52')
    onames <- c('lM','l25','l45','l65','l75','l85')
    ierrnames <- c('errmgspk','err74','err64',
                   'err76','err65','err86','err52')
    lsdat <- unlist(log(samples[i,inames]))
    errlsdat <- samples[i,ierrnames] # relative errors
    spk <- samples[i,'spk']
    ispk <- which(spk%in%spikes[,'name'])
    lM <- lsdat['mgspk'] + ifelse(conc,log(spikes[ispk,'pmg205']),0)
    lr25 <- -lsdat['r52']
    lr45 <- lsdat['r65']-lsdat['r64']
    lr65 <- lsdat['r65']
    lr75 <- lsdat['r76']+lsdat['r65']
    lr85 <- lsdat['r86']+lsdat['r65']
    lsmp <- c(lM,lr25,lr45,lr65,lr75,lr85)
    names(lsmp) <- onames
    E <- diag(errlsdat)^2
    J <- matrix(0,nrow=6,ncol=7)
    colnames(J) <- inames
    rownames(J) <- onames
    J['lM','mgspk'] <- 1
    J['l25','r52'] <- -1
    J['l45','r65'] <- 1
    J['l45','r64'] <- -1
    J['l65','r65'] <- 1
    J['l75','r76'] <- 1
    J['l75','r65'] <- 1
    J['l85','r86'] <- 1
    J['l85','r65'] <- 1
    covlsmp <- J %*% E %*% t(J)
    spk <- samples[i,'spk']
    lspk <- log(spikes[ispk,'r52'])
    # placeholder for spikes[spk[1],'err52']
    errlspk <- 0
    list(lsmp=lsmp,covlsmp=covlsmp,lspk=lspk,errlspk=errlspk)
}
getsample <- function(i,samples,spikes,conc=TRUE){
    a <- getaliquot(i=i,samples=samples,spikes=spikes,conc=conc)
    spk <- samples[i,'spk']
    ispk <- which(spk%in%spikes[,'name'])
    out <- lxstar(lM=a$lsmp['lM'], 
                  lr45=a$lsmp['l45'],
                  lr65=a$lsmp['l65'],
                  lr75=a$lsmp['l75'],
                  lr85=a$lsmp['l85'],
                  lr25=a$lsmp['l25'],
                  lr25t=-log(spikes[ispk,'r52']))
    out
}

getE <- function(i=1,samples,ablk,spikes){
    lsmp <- getaliquot(i,samples,spikes)
    Es <- lsmp$covlsmp[c('lM','l25','l45','l65','l75'),
                       c('lM','l25','l45','l65','l75')]
    Eb <- ablk$covlblk[c('l4','l6','l7'),c('l4','l6','l7')]
    E <- matrix(0,9,9)
    E[1:3,1:3] <- Eb
    E[4:8,4:8] <- Es
    if (lsmp$errlspk != ablk$errlspk){
        stop("The spike and sample don't use the same blank.")
    } else {
        E[9,9] <- lsmp$errlspk^2
    }
    J <- matrix(0,nrow=6,ncol=9)
    J[1:3,1:3] <- diag(3)
    J[4:6,4] <- 1
    J[4:6,5] <- -beta(c(204,206,207))
    J[4:6,6:8] <- diag(3)
    J[1:3,9] <- beta(c(204,206,207))
    J[4:6,9] <- J[1:3,9]
    J %*% E %*% t(J)
}

LL <- function(p,i=1,lsmp,lblk,E){
    c4b <- p[1]
    c6b <- p[2]
    c7b <- p[3]
    c6 <- p[4]
    c46 <- p[5]
    c76 <- p[6]
    D <- rep(NA,6)
    D[1] <- lblk['l4'] - c4b
    D[2] <- lblk['l6'] - c6b
    D[3] <- lblk['l7'] - c7b
    D[4] <- lsmp[1,'l4'] - log(exp(c46+c6) + exp(c4b))
    D[5] <- lsmp[1,'l6'] - log(exp(c6b) + exp(c6))
    D[6] <- lsmp[1,'l7'] - log(exp(c76+c6) + exp(c7b))
    stats::mahalanobis(D,center=FALSE,cov=E)/2
}

init <- function(i,lsmp,lblk){
    c4b <- lblk['l4']
    c6b <- lblk['l6']
    c7b <- lblk['l7']
    c6 <- lsmp[1,'l6']
    c46 <- lsmp[1,'l4'] - lsmp[1,'l6']
    c76 <- lsmp[1,'l7'] - lsmp[1,'l6']
    out <- c(c4b,c6b,c7b,c6,c46,c76)
    names(out) <- c('l4b','l6b','l7b','l6','l46','l76')
    out
}

#' @title Pb-Pb data processing
#' @description High level functionthat takes samples, blanks and
#'     spikes as input and produces a table of Pb/Pb ratios,
#'     uncertainties and error correlations as output
#' @param samples data frame with sample data. For an example, see
#'     \code{system.file("samples1.csv",package="PbPbRedux")}
#' @param blanks data frame with blank data. For an example, see
#'     \code{system.file("blanks1.csv",package="PbPbRedux")}
#' @param spikes data frame with spike data. For an example, see
#'     \code{system.file("spikes.csv",package="PbPbRedux")}
#' @param cblanks data frame with replicate blank data. For an
#'     example, see
#'     \code{system.file("blanks2.csv",package="PbPbRedux")}
#' @param ierr indicates whether the analytical uncertainties are
#'     reported as: 
#' 
#' \code{1}: 1\eqn{\sigma} absolute uncertainties.
#' 
#' \code{2}: 2\eqn{\sigma} absolute uncertainties.
#' 
#' \code{3}: 1\eqn{\sigma} relative uncertainties (\eqn{\%}).
#' 
#' \code{4}: 2\eqn{\sigma} relative uncertainties (\eqn{\%}).
#' 
#' @return a table with Pb concentrations, Pb/Pb ratios and error
#'     correlations
#' @examples
#' library(PbPbRedux)
#' 
#' s1 <- system.file("samples1.csv",package="PbPbRedux")
#' s2 <- system.file("samples2.csv",package="PbPbRedux")
#' b1 <- system.file("blanks1.csv",package="PbPbRedux")
#' b2 <- system.file("blanks2.csv",package="PbPbRedux")
#' sp <- system.file("spikes.csv",package="PbPbRedux")
#'
#' spikes <- read.csv(sp,header=TRUE)
#' 
#' # example 1: all samples use the same blank:
#' samples <- read.csv(s1,header=TRUE)
#' blanks <- read.csv(b1,header=TRUE)
#' tab <- process(samples,blanks,spikes)
#'
#' # example 2: each aliquot has its own blank:
#' samples <- read.csv(s2,header=TRUE)
#' blanks <- read.csv(b2,header=TRUE)
#' tab <- process(samples,blanks,spikes)
#'
#' # example 3: individual blanks with shared covariance matrix:
#' samples <- read.csv(s2,header=TRUE)
#' cblanks <- read.csv(b1,header=TRUE)
#' blanks <- read.csv(b2,header=TRUE)
#' tab <- process(samples,blanks,spikes,cblanks)
#' 
#' @export
process <- function(samples,blanks,spikes,cblanks,ierr=4){
    cb <- !missing(cblanks)
    SM <- ierradj(samples,2*(2:9),ierr=ierr)
    BL <- ierradj(blanks,2*(1:4)+1,ierr=ierr)
    SP <- spikes
    ns <- nrow(samples)
    cnames <- c('4/6','err[4/6]','7/6','err[7/6]','rho','pgPb','pgPb(blk)')
    if (cb){ # each aliquot has its own blank but covariance matrix is shared
        CB <- ierradj(cblanks,2*(1:4)+1,ierr=ierr)
        cblk <- avgblank(i=1:ns,blanks=CB,spikes=SP)$covlblk
    }
    out <- matrix(NA,nrow=ns,ncol=length(cnames))
    colnames(out) <- cnames
    rownames(out) <- SM$Label
    for (ii in 1:nrow(SM)){
        ablk <- avgblank(blanks=BL,blk=SM[ii,'blk'],spikes=SP,spk=SM[ii,'spk'])
        if (cb) ablk$covlblk <- cblk
        lsmp <- getsample(i=ii,SM,SP)
        E <- getE(ii,SM,ablk,SP)
        pinit <- init(i=ii,lsmp=lsmp,lblk=ablk$lblk)
        fit <- optim(pinit,fn=LL,method='BFGS',i=ii,lsmp=lsmp,
                     lblk=ablk$lblk,E=E,hessian=TRUE)
        J <- diag(exp(fit$par))
        H <- IsoplotR:::nearPD(fit$hessian)
        covmat <- J %*% solve(H) %*% t(J)
        cormat <- cov2cor(covmat)
        out[ii,'pgPb'] <- sum(exp(lsmp)*c(204,206:208))
        out[ii,'pgPb(blk)'] <- sum(exp(ablk$lblk)*c(204,206:208))
        out[ii,c('4/6','7/6')] <- exp(fit$par[5:6])
        out[ii,'err[4/6]'] <- sqrt(covmat[5,5])
        out[ii,'err[7/6]'] <- sqrt(covmat[6,6])
        out[ii,'rho'] <- cormat[5,6]
    }
    out <- oerradj(out,c(1,3),oerr=ierr)
    out
}
