#IsotopeR model. Bayesian Inference on stable istope analysis.
#code based on Semmens et al 2009
#see IsotopeR instructions.doc for details on running the model.
#see Hopkins and Ferguson 2011 for details on the model structure.
#Jake Ferguson: last updated 4/25/12


IsoWrapper <- function(Mixtures="Necessary File", Sources="Necessary File", Concentrations="Optional File", Measurement.Error="Optional File", Discrimination.Error="Optional File", Digestibility.Factor="Optional File", output.name="SampleOutput.Rdata", mcmc.chains=3, mcmc.burn=1e3, mcmc.chainLength=1e3, mcmc.thin=10, plot.mixing.observations=TRUE, plot.mixing.estimates=TRUE, plot.dietary.source.contributions=TRUE) {
    require(runjags)
#   source('IsotopeRModelsNoGroups.R')
    
    mcmc.chainLength    <- as.integer(mcmc.chainLength+mcmc.burn) #total number of iterations per chain (includes burnin)

    #name and location of the model file passed to JAGS
    model.loc   <- "IsotopeR.txt" 
    #parameters to be returned by JAGS
    jags.params <- c("mu.source", "sd.source", "rho.source", "mu.conc", "sd.conc", "mu.mix", "p", "p.pop", "sd.me","sd.res")
               
    file.flag <- NULL
    noconc.flag = 0
    nome.flag = 0
    nodiscrim.flag = 0
    nodigest.flag=0
    
    #reads in the files
    X           <- try(as.matrix(read.table(Mixtures, sep='\t', header=TRUE)), silent=TRUE) #mixture data file
    if(class(X) == 'try-error') {stop("Mixture file not found")}
    if(dim(X)[2] == 1) {
        X           <- as.matrix(read.table(Mixtures, sep=',', header=TRUE)) #mixture data file
    }
    
    sources     <- try((read.table(Sources,sep='\t',header=TRUE)), silent=TRUE) #source data
    if(class(sources) == 'try-error') {stop("Source file not found")}
    if(dim(sources)[2] == 1) { 
        sources     <- (read.table(Sources,sep=',',header=TRUE)) #source data
    }
    
    D   <- NA
    if(Concentrations == 'Optional File') { 
        noconc.flag = 1
        file.flag       <- paste(file.flag,"noconc",sep='') 
        jags.rem    <- which(jags.params == 'sd.conc')
        jags.params <- jags.params[-jags.rem]
    } else {
        D           <- try((read.table(Concentrations,sep='\t',header=TRUE)), silent=TRUE) #source concentration data
        if(dim(D)[2] == 1) {
            D           <- try(read.table(Concentrations,sep=',',header=TRUE), silent=TRUE) #source concentration data
        }
    }
    
    Z <- NA
    if(Measurement.Error == 'Optional File') { 
        nome.flag = 1
        file.flag   <- paste(file.flag,"nome",sep='') 
        jags.rem    <- which(jags.params == 'sd.me')
        jags.params <- jags.params[-jags.rem]
    } else {
        Z           <- try(as.matrix(read.table(Measurement.Error,sep='\t',header=TRUE)), silent=TRUE) #file for measurement error
        if(class(Z) == 'try-error') {stop("Measurement Error file not found")}
        if(dim(Z)[2] == 1) { 
            Z           <- as.matrix(read.table(Measurement.Error, sep=',', header=TRUE)) #file for measurement error
        }
    }
    
    discrim.sd <- NA
    if(Discrimination.Error == 'Optional File') {
        nodiscrim.flag = 1
        file.flag   <- paste(file.flag,"nodiscrim",sep='')
    } else {
        discrim.sd      <- try((read.table(Discrimination.Error,sep='\t',header=TRUE)), silent=TRUE) #file with the standard deviation of discrimination in sources
        if(class(discrim.sd) == 'try-error') {stop("Discrimination Error file not found")}
        if(dim(discrim.sd)[2] == 1) {
            discrim.sd      <- (read.table(Discrimination.Error,sep=',',header=TRUE)) #file with the standard deviation of discrimination in sources
        }
    }

    digest <- NA
    if(Digestibility.Factor == 'Optional File') {
        nodigest.flag = 1
#         file.flag   <- paste(file.flag,"nodigest",sep='')
    } else {
        digest      <- try((read.table(Digestibility.Factor,sep='\t',header=TRUE)), silent=TRUE) #file with the standard deviation of discrimination in sources
        if(class(digest) == 'try-error') {stop("Digestibility factor file not found")}
        if(dim(digest)[2] == 1) {
            digest      <- (read.table(Digestibility.Factor,sep=',',header=TRUE)) #file with digestion
        }
    }
    
    
    if(file.flag == "") { curr.model <- IsotopeR}
    if(file.flag == "noconcnomenodiscrim") { curr.model <- IsotopeRnoconcnomenodiscrim}
    if(file.flag == "noconc") { curr.model <- IsotopeRnoconc}
    if(file.flag == "noconcnodiscrim") { curr.model <- IsotopeRnoconcnodiscrim}
    if(file.flag == "noconcnome") { curr.model <- IsotopeRnoconcnome}
    if(file.flag == "nodiscrim") { curr.model <- IsotopeRnodiscrim}
    if(file.flag == "nome") { curr.model <- IsotopeRnome}
    if(file.flag == "nomenodiscrim") { curr.model <- IsotopeRnomenodiscrim}

    ######Do not mess below this unless you understand what will happen#######
    #extract useful info from the read in data to pass to JAGS
    N       <- dim(X)[1] #number of individuals in the sample
    num.iso     <- dim(X)[2]-2 #number of isotopes in the sample
    num.sources <- nlevels(as.factor(sources[,num.iso+1]))

    #prior diet proportions
    alpha <- rep(1,num.sources)/num.sources

    #center Z (observation error) around 0
    Nz <- NA
    if(!nome.flag) {
        Z[,1] <- Z[,1] - mean(Z[,1])
        Z[,2] <- Z[,2] - mean(Z[,2])
        Nz <- dim(Z)[1]
    }

    
    #prior parameters for concentrations
    dmu.prior.tau <- diag(num.iso)*1/1000
    dmu.prior.mu <- matrix(c(50,1),num.iso,num.sources)

    #population prior
    alpha.clr <- log(alpha/prod(alpha)^(1/length(alpha))) #transform onto CLR scale

    #priors for concentration covariance matrix
    D.R <- diag(num.iso)*1/1000

    #measurement error covariance matrix and mean
    tauZ    <- diag(num.iso)
    muz     <- rep(0,num.iso)
    
    #puts sources into an array- makes it easier to use in JAGS (this is just processing stuff and not very important to understand)   
    source.ss <- vector('numeric', num.sources)
    counter <- 1

    if(class(sources$source) == 'factor') {
        for(i in levels(sources$source)) {    
            source.ss[counter] <- length(which(sources$source == i))
            counter <- counter+1        
        }    
        names(source.ss) <- levels(sources$source)
    
        source.array <- array(NA, c(max(source.ss), num.iso,num.sources))
        counter <- 1

        for(i in levels(as.factor(sources$source))) {
            source.array[1:source.ss[counter],,counter] <- as.matrix(sources[which(sources$source == i),1:num.iso])
            counter <- counter+1        
        }
    } else {
        for(i in unique(sources$source)) {    
            source.ss[counter] <- length(which(sources$source == i))
            counter <- counter+1        
        }    
        names(source.ss) <- unique(sources$source)
        source.array <- array(NA, c(max(source.ss), num.iso,num.sources))
        counter <- 1

       for(i in unique(sources$source)) {            
            source.array[1:source.ss[counter],,counter] <- as.matrix(sources[which(sources$source == i),1:num.iso])
            counter <- counter+1        
        }
    }
    
    #read in subgroups
    num.subgroups <- sources[,num.iso + 2]
    if(levels(as.factor(num.subgroups)) != 1) {stop('multiple groups not implemented in this version')}

    
    #prior parameters for sources  
    mu.prior.mu <- c(0, 0) #apply(sources[,1:num.iso], 2, mean)
    mu.prior.cov <- solve(diag(num.iso)*100) #source mean prior covariance matrix
    
    cd.array <- NA
    cd.ss <- NA
    if(!noconc.flag) {
        #bu ild array of isotope concentrations and apply digestability
        cd.ss <- vector('numeric', num.sources)
        counter <- 1
        if(class(D$source) =='factor') {
            for(i in levels(D$source)) {
                cd.ss[counter] <- length(which(D$source == i))
                counter <- counter+1
            }
            cd.array <- array(NA, c(max(cd.ss), num.iso, num.sources))
            counter <- 1
            for(i in levels(D$source)) {
                cd.array[1:cd.ss[counter],,counter] <- as.matrix(D[which(D$source == i),1:num.iso])
                if(!nodigest.flag) {cd.array[,,counter] <- t(t(cd.array[,,counter])*digest[counter,1:num.iso]) }
                counter <- counter+1
            }
        } else {
           for(i in unique(D$source)) {
                cd.ss[counter] <- length(which(D$source == i))
                counter <- counter+1
            }
            cd.array <- array(NA, c(max(cd.ss), num.iso, num.sources))
            counter <- 1
            for(i in unique(D$source)) {
                cd.array[1:cd.ss[counter],,counter] <- as.matrix(D[which(D$source == i),1:num.iso])
                if(!nodigest.flag) {cd.array[,,counter] <- t(t(cd.array[,,counter])*digest[counter,1:num.iso]) }
                counter <- counter+1
            }
        }
    }
    #builds discrimination variation matrices (added to estimators of the standard deviation in order to account for discrimination variation)
    discrimvar.mat <- NA
    if(!nodiscrim.flag) {
        discrim.var  <- (apply(discrim.sd[,-(num.iso+1)]^2,2, sum) )/dim(discrim.sd)[1]^2
        discrimvar.mat  <- diag(discrim.var)
    }
    #gets number of individual observations
    num.inds <- nlevels(as.factor(X[,4]))
    ind.counts <- vector('numeric',num.inds)
    ind.counts.index    <- matrix(NA, num.inds, dim(X)[1]-num.inds+1)
    
    for(i in 1:num.inds) {
        ind.counts[i] <- length(which(i == X[,4]))
        ind.counts.index[i,] <- which(X[,4] == i)
    }

    #individual observation matrix
    ind.array <- array(NA,c(num.iso,num.inds,max(ind.counts)))
    for(i in 1:num.inds) {
        ind.array[1:num.iso,i,] <- t(X[which(X[,4]==i),1:num.iso])
    }
    #individual observation id's
    
    N <- num.inds
    jags.data = list( "muz", "ind.counts", "ind.array", "ind.counts.index", "N", "num.sources", "num.iso", "Z", "dmu.prior.mu", "Nz", "source.array", "cd.array", "mu.prior.mu", "mu.prior.cov", "dmu.prior.tau", "alpha.clr", "source.ss", "cd.ss", "discrimvar.mat")  #data passed into JAGS
   
   jags.dump <- list(muz=muz, ind.counts=ind.counts, ind.array=ind.array, ind.counts.index=ind.counts.index, N=N, num.sources=num.sources, num.iso=num.iso, Z=Z, dmu.prior.mu=dmu.prior.mu, Nz=Nz, source.array=source.array, cd.array=cd.array, mu.prior.mu=mu.prior.mu, mu.prior.cov=mu.prior.cov, dmu.prior.tau=dmu.prior.tau, alpha.clr=alpha.clr, source.ss=source.ss, cd.ss=cd.ss, discrimvar.mat=discrimvar.mat)

    if(noconc.flag) {
        jags.rem    <- which(jags.data  == 'cd.array')
        jags.data <- jags.data[-jags.rem]
        jags.rem    <- which(jags.data  == 'cd.ss')
        jags.data <- jags.data[-jags.rem]
    }
    if(nome.flag) {
        jags.rem    <- which(jags.data  == 'Z' | jags.data == 'muZ' | jags.data == 'tauZ' | jags.data== 'Nz')
        jags.data <- jags.data[-jags.rem]
    }
    if(nodiscrim.flag) {
        jags.rem    <- which(jags.data  == 'discrimvar.mat')
        jags.data <- jags.data[-jags.rem]
    }
    #function used to initialize parameters
    jags.inits <- list(dmu.prior.mu=dmu.prior.mu, mu.prior.mu=mu.prior.mu, p.transform=runif(num.sources), region.sig=0.5, ind.sig=0.5, p.ind = matrix(runif(N*num.sources), N, num.sources))
write.table(jags.data, file='temp.txt')
    jags.out <- run.jags(model=curr.model, monitor=jags.params, data=jags.dump, n.chains=mcmc.chains, burnin=mcmc.burn, sample= mcmc.chainLength, thin=mcmc.thin, check.conv=TRUE, plots=FALSE, monitor.deviance=FALSE, silent.jags=TRUE)
write.table(jags.data, file='temp2.txt')
    r.est <- jags.out$psrf$psrf[,1]
    jags.output.mat <- cbind(jags.out$summary$statistics[,1:2], jags.out$summary$quantiles, r.est)
    
    save(jags.out, X, file=output.name) #save image file with all the info from the MCMC 
    write.table(jags.output.mat, file=paste(strsplit(output.name, ".Rdata")[[1]],'.txt',sep=''))

    
    if(plot.mixing.observations) {
#        source("Plot_jags.r")
        X11()
        if(num.sources == 3) { Tri.plots(jags.out,X,plot.ind=TRUE, me.flag=!nome.flag) } else {
        if(num.sources == 2) { Bi.plots(jags.out,X,plot.ind=TRUE, me.flag=!nome.flag) } else {stop("No plotting method available for this many isotopes")}
        }
    }

    if(plot.mixing.estimates) {
        X11()
        if(num.sources == 3) { Tri.plots(jags.out, X, plot.mix=TRUE, me.flag=!nome.flag) } else {
        if(num.sources == 2) { Bi.plots(jags.out, X,  plot.mix=TRUE, me.flag=!nome.flag) } else {stop("No plotting method available for this many isotopes")}
        }    
    }
    
    if(plot.dietary.source.contributions) {
        X11()
        curves.plot(jags.1=jags.out, num.sources=num.sources, num.chains=mcmc.chains, color=TRUE, individuals=N,  xlab.vec=unique(sources[,num.iso+1]))
    }

    print(jags.output.mat)
    return(jags.out)

}

IsotopeR    <- function() {
    options(warn=-1) 
require(fgui) #used to run JAGS from R

output <- guiv(IsoWrapper, argFilter = list(Mixtures="{{} {.csv}}", Sources="{{} {.csv}}", Concentrations="{{} {.csv}}", Measurement.Error="{{} {.csv}}",   Discrimination.Error="{{} {.csv}}", Digestibility.Factor="{{} {.csv}}", output.name="{{} {.Rdata}}"), argText = list(mcmc.chains="number of chains", mcmc.burn="MCMC burnin", mcmc.chainLength="MCMC runs", mcmc.thin="thinning rate"), argOption = list(plot.mixing.observations=c("TRUE","FALSE"), plot.dietary.source.contributions= c("TRUE", "FALSE"), plot.mixing.estimates=c('TRUE','FALSE')), closeOnExec=TRUE, title="IsotopeR v0.1", exec="Run IsotopeR", helpsFunc="IsoWrapper")
     options(warn=0)
}
