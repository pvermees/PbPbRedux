beta <- function(mx){
    (log(205)-log(mx))/(log(205)-log(202))
}

lxstar <- function(lM,lr45,lr65,lr75,lr25,lr25t){
    if (missing(lr25t)) fract <- -lr25
    else fract <- lr25t-lr25
    l4 <- lr45 + lM + beta(204)*fract
    l6 <- lr65 + lM + beta(206)*fract
    l7 <- lr75 + lM + beta(207)*fract
    out <- cbind(l4,l6,l7)
    colnames(out) <- c('l4','l6','l7')
    out
}

avgblank <- function(blanks,spikes,spk){
    i <- which(blanks[,'spk']%in%spk)
    lbdat <- log(blanks[i,c('mgspk','r52','r54','r56','r57'),drop=FALSE])
    lblk <- lxstar(lM=lbdat[,'mgspk'] +
                       log(spikes[spk,'pmg205']), # pmg205 optional
                   lr45=-lbdat[,'r54'], # flip sign
                   lr65=-lbdat[,'r56'], # flip sign
                   lr75=-lbdat[,'r57'], # flip sign
                   lr25=-lbdat[,'r52'], # flip sign
                   lr25t=-log(spikes[spk,'r52'])) # flip sign
    lblkprime <- lxstar(lM=lbdat[,'mgspk'] +
                            log(spikes[spk,'pmg205']), # pmg205 optional
                        lr45=-lbdat[,'r54'], # flip sign
                        lr65=-lbdat[,'r56'], # flip sign
                        lr75=-lbdat[,'r57'], # flip sign
                        lr25=-lbdat[,'r52']) # flip sign
    mlblk <- colMeans(lblkprime)
    if (nrow(lblkprime)<2){
        covlblk <- matrix(0,3,3)
    } else {
        covlblk <- cov(lblkprime)
    }
    lspk <- -log(spikes[spk,'r52'])[1]
    errlspk <- 0 # placeholder for spikes[spk[1],'err52']/200
    list(lblk=mlblk,covlblk=covlblk,lspk=lspk,errlspk=errlspk)
}

# get x/5 data from a particular aliquot:
getaliquot <- function(i,samples,spikes){
    inames <- c('mgspk','r74','r64','r76','r65','r52')
    onames <- c('lM','l25','l45','l65','l75')
    lsdat <- unlist(log(samples[i,inames]))
    errlsdat <- samples[i,c('errmgspk','err74','err64',
                            'err76','err65','err52')]/200
    spk <- samples[i,'spk']
    lM <- lsdat['mgspk'] + log(spikes[spk,'pmg205']) # pmg205 optional
    lr25 <- -lsdat['r52']
    lr45 <- lsdat['r65']-lsdat['r64']
    lr65 <- lsdat['r65']
    lr75 <- lsdat['r76']+lsdat['r65']
    lsmp <- c(lM,lr25,lr45,lr65,lr75)
    names(lsmp) <- onames
    E <- diag(errlsdat)^2
    J <- matrix(0,nrow=5,ncol=6)
    colnames(J) <- inames
    rownames(J) <- onames
    J['lM','mgspk'] <- 1
    J['l25','r52'] <- -1
    J['l45','r65'] <- 1
    J['l45','r64'] <- -1
    J['l65','r65'] <- 1
    J['l75','r76'] <- 1
    J['l75','r65'] <- 1
    covlsmp <- J %*% E %*% t(J)
    spk <- samples[i,'spk']
    lspk <- log(spikes[spk,'r52'])
    # placeholder for spikes[spk[1],'err52']/200
    errlspk <- 0
    list(lsmp=lsmp,covlsmp=covlsmp,lspk=lspk,errlspk=errlspk)
}
getsample <- function(i,samples,spikes){
    a <- getaliquot(i=i,samples=samples,spikes=spikes)
    spk <- samples[i,'spk']
    out <- lxstar(lM=a$lsmp['lM'], # pmg205 optional
                  lr45=a$lsmp['l45'],
                  lr65=a$lsmp['l65'],
                  lr75=a$lsmp['l75'],
                  lr25=a$lsmp['l25'],
                  lr25t=-log(spikes[spk,'r52']))
    out
}

getE <- function(i=1,samples,ablk,spikes){
    lsmp <- getaliquot(i,samples,spikes)
    Es <- lsmp$covlsmp # lx5 data
    Eb <- ablk$covlblk # lx data
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
    (D %*% solve(E) %*% D)
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

process <- function(samples,blanks,spikes){
    ns <- nrow(samples)
    out <- matrix(NA,nrow=ns,ncol=5)
    colnames(out) <- c('4/6','s[4/6]','7/6','s[7/6]','rho')
    rownames(out) <- samples$Label
    for (i in 1:nrow(samples)){
        print(i)
        ablk <- avgblank(blanks,spikes,spk=samples[i,'spk']) # with covariance matrix
        lsmp <- getsample(i=i,samples,spikes) # without covariance matrix
        E <- getE(i,samples,ablk,spikes)
        pinit <- init(i=i,lsmp=lsmp,lblk=ablk$lblk)
        fit <- optim(pinit,fn=LL,i=i,lsmp=lsmp,
                     lblk=ablk$lblk,E=E,hessian=TRUE)
        J <- diag(exp(fit$par))
        covmat <- J %*% solve(fit$hessian) %*% t(J)
        cormat <- cov2cor(covmat)
        out[i,c('4/6','7/6')] <- exp(fit$par[5:6])
        out[i,'s[4/6]'] <- sqrt(covmat[5,5])
        out[i,'s[7/6]'] <- sqrt(covmat[6,6])
        out[i,'rho'] <- cormat[5,6]
    }
    out
}
