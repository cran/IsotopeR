#IsotopeR model. Bayesian Inference on stable istope analysis.
#code based on Semmens et al 2009
#see IsotopeR instructions.doc for details on running the model.
#see Hopkins and Ferguson 2011 for details on the model structure.
#Jake Ferguson: last updated 7/1/11

IsoWrapper <- function(Mixtures="Necessary File", Sources="Necessary File", Concentrations="Optional File", Measurement.Error="Optional File", Discrimination.Error="Optional File", Digestibility.Factor="Optional File", output.name="SampleOutput.Rdata", mcmc.chains=3, mcmc.burn=1e3, mcmc.chainLength=1e3, mcmc.thin=100, plot.observations=TRUE,  plot.mixing.estimates=TRUE, plot.dietary.source.contributions=TRUE, color.plots=TRUE, run.parallel = TRUE) {
    require(runjags)
    
    mcmc.chainLength    <- as.integer(mcmc.chainLength+mcmc.burn) #total number of iterations per chain (includes burnin)

    #name and location of the model file passed to JAGS
    model.loc   <- "IsotopeR.txt" 

    #parameters to be returned by JAGS
    jags.params <- c("mu.source", "sd.source", "rho.mat", "mu.conc", "sd.conc", "mu.mix", "p", "p.pop", "sd.me","sd.res")
    file.flag <- ""
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
    cd.mat <- NA
    subcd.vec <- NA
    subcd.samples <- NA
	if(Concentrations == 'Optional File') { 
        noconc.flag = 1
        file.flag       <- paste(file.flag,"noconc",sep='') 
        
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

    options(warn=-1)
    #extract useful info from the read in data to pass to JAGS
    N       				<- dim(X)[1] #number of individuals in the sample
    num.iso     		<- dim(X)[2]-2 #number of isotopes in the sample
    num.sources 	<- nlevels(as.factor(sources[,num.iso+1]))

	num.groups 		<- nlevels(as.factor(X[,num.iso+1]))
	groupnum.mat 	<- matrix(NA,num.groups,2)
	for(i in 1:num.groups) {
		groupnum.mat[i,1] <- min(which(X[,num.iso+1]==i))
		groupnum.mat[i,2] <- max(which(X[,num.iso+1]==i))
	}
	#determine the proper model to run

	if(file.flag == "") { 
		if(num.groups ==1) { curr.model <- IsotopeRfull} else {curr.model <- IsotopeRfullgroup}
	}
    if(file.flag == "noconcnomenodiscrim") { 
		if(num.groups ==1) { curr.model <- IsotopeRnoconcnomenodiscrim} else { curr.model <- IsotopeRnoconcnomenodiscrimgroup}
	}
    if(file.flag == "noconc") { 
		if(num.groups ==1) { curr.model <- IsotopeRnoconc } else { curr.model <- IsotopeRnoconcnomenodiscrimgroup}
	}
    if(file.flag == "noconcnodiscrim") { 
		if(num.groups ==1) { curr.model <- IsotopeRnoconcnodiscrim} else {curr.model <- IsotopeRnoconcnodiscrimgroup}
	}
    if(file.flag == "noconcnome") { 
		if(num.groups == 1) {curr.model <- IsotopeRnoconcnome} else {curr.model <- IsotopeRnoconcnomegroup}
	}
    if(file.flag == "nodiscrim") { 
		if(num.groups == 1) { curr.model <- IsotopeRnodiscrim } else {curr.model <- IsotopeRnodiscrimgroup}
	}
    if(file.flag == "nome") { 
		if(num.groups == 1) { curr.model <- IsotopeRnome } else {curr.model <- IsotopeRnomegroup }
	}
    if(file.flag == "nomenodiscrim") { 
		if(num.groups == 1) { curr.model <- IsotopeRnomenodiscrim } else {curr.model <- IsotopeRnomenodiscrimgroup }
	}

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
    dmu.prior.mu <- matrix(100,num.iso,num.sources)

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

        for(i in levels(sources[,num.iso+2])) {    
            source.ss[counter] <- length(which(sources[,num.iso+2] == i))
            counter <- counter+1        
        }    
        names(source.ss) <- levels(sources$source)
    
		##get array indices for the different subsources
		subsources <- sources[,num.iso+2]
		subsource.vec <- nlevels(as.factor(sources[,num.iso+1]))
		counter <- 1
		for(i in levels(as.factor(sources[,num.iso+1]))) {
			subsource.vec[counter] <- nlevels(as.factor(sources[which(sources[,num.iso+1] == i), num.iso+2]))
			counter <- counter+1
		}
		subsource.samples <-array(NA, c(num.sources, max(subsource.vec), 2))
		source.counter <- 1		
		for(i in levels(as.factor(sources[,num.iso+1]))) {
 			curr.source.index <- which(sources[,num.iso+1] == i)
			subsource.counter <- 1
			for(j in levels(as.factor(sources[curr.source.index, num.iso+2]))) {
				curr.subsource.index <- which(sources[curr.source.index, num.iso+2] == j)
				subsource.samples[source.counter, subsource.counter, ] <- 	curr.source.index[1] + curr.subsource.index[c(1, length(curr.subsource.index))] - 1				
				
				subsource.counter <- subsource.counter +1
				
			}
			source.counter <- source.counter + 1
		}
        counter <- 1
        source.mat <- as.matrix(sources[,1:num.iso])
    
    #prior parameters for sources  
    mu.prior.mu <- rep(0,num.iso)#c(0, 0) #apply(sources[,1:num.iso], 2, mean)
    mu.prior.cov <- solve(diag(num.iso)*100) #source mean prior covariance matrix
    cd.array <- NA

    if(!noconc.flag) {
        #bu ild array of isotope concentrations and apply digestability
#         if(class(D[,num.iso+1]) =='factor') {
        cd.mat <- as.matrix(D[,1:num.iso])

         ##get array indices for the different subsource concentrations
		subcd <- D[,num.iso+2]
		subcd.vec <- nlevels(as.factor(D[,num.iso +1]))
		counter <- 1
		for(i in levels(as.factor(D[,num.iso+1]))) {
			subcd.vec[counter] <- nlevels(as.factor(D[which(D[,num.iso+1] == i), num.iso+2]))
			counter <- counter+1
		}

		subcd.samples <-array(NA, c(num.sources, max(subcd.vec), 2))
		source.counter <- 1
		for(i in levels(as.factor(D[,num.iso+1]))) {
 			curr.source.index <- which(D[,num.iso+1] == i)

			subsource.counter <- 1
			for(j in levels(as.factor(D[curr.source.index, num.iso+2]))) {
				curr.subsource.index <- which(D[curr.source.index, num.iso+2] == j)
				
				subcd.samples[source.counter, subsource.counter, ] <- 	curr.source.index[1] + curr.subsource.index[c(1, length(curr.subsource.index))] - 1
				subsource.counter <- subsource.counter +1
			}
			source.counter <- source.counter + 1
		}
    }

	#builds discrimination variation matrices (added to estimators of the standard deviation in order to account for discrimination variation)
    discrimvar.mat <- NA

    if(!nodiscrim.flag) {
        discrim.var  <- (apply(discrim.sd[,-(num.iso+1)]^2,2, sum) )/dim(discrim.sd)[1]^2
        discrimvar.mat  <- diag(discrim.var)
    }
    #gets number of individual observations
	dim.x	<- dim(X)
    num.inds <- nlevels(as.factor(X[,dim.x[2]]))
    ind.counts <- vector('numeric',num.inds)
    for(i in 1:num.inds) {
        ind.counts[i] <- length(which(i == X[,dim.x[2]]))
    }
    #individual observation matrix

    ind.array <- array(NA,c(num.iso,num.inds,max(ind.counts)))
    for(i in 1:num.inds) {
        ind.array[1:num.iso,i,] <- t(X[which(X[,dim.x[2]]==i), 1:num.iso])
    }

# 	print(Tau.mat)
# 	tau.mat	<- diag(num.iso)
# 	tau.df		<- num.iso+100
	rho.flag	<- ifelse(num.iso ==2, 1, 0)
	#individual observation id's
    N <- num.inds
   jags.dump <- list(muz=muz, ind.counts=ind.counts, ind.array=ind.array,N=N, num.sources=num.sources, num.iso=num.iso, Z=Z, dmu.prior.mu=dmu.prior.mu, Nz=Nz, mu.prior.mu=mu.prior.mu, mu.prior.cov=mu.prior.cov, dmu.prior.tau=dmu.prior.tau, alpha.clr=alpha.clr,  discrimvar.mat=discrimvar.mat, subsource.vec=subsource.vec, subsource.samples=subsource.samples, source.mat=source.mat, cd.mat = cd.mat, subcd.vec=subcd.vec, subcd.samples=subcd.samples, num.groups=num.groups, groupnum.mat=groupnum.mat, rho.flag=rho.flag)

   if(noconc.flag) {
        jags.rem    <- which( names(jags.dump)  == 'cd.mat' |  names(jags.dump)  == 'dmu.prior.mu' |  names(jags.dump)  == 'dmu.prior.tau' |  names(jags.dump)  == 'subcd.samples' |  names(jags.dump)  == 'subcd.vec')
        jags.dump <- jags.dump[-jags.rem]
        
       jags.rem <- which(jags.params == "mu.conc" | jags.params=="sd.conc") 
	   jags.params <- jags.params[-jags.rem]
	
    }
    if(nome.flag) {
        jags.rem    <- which(names(jags.dump)  == 'Z' | names(jags.dump) == 'muz' | names(jags.dump) == 'tauZ' | names(jags.dump)== 'Nz')
        jags.dump <- jags.dump[-jags.rem]
        
       jags.rem <- which(jags.params == "sd.me") 
	   jags.params <- jags.params[-jags.rem]
    }
    if(nodiscrim.flag) {
        jags.rem    <- which(names(jags.dump)  == 'discrimvar.mat')
        jags.dump <- jags.dump[-jags.rem]
    }    
    if(num.groups <= 1) {
		jags.rem    <- which(names(jags.dump)  == 'groupnum.mat' | names(jags.dump)  == 'num.groups')
        jags.dump <- jags.dump[-jags.rem]        
	} else {
		jags.params <- c(jags.params,"p.group") 		
	}
	if(num.iso != 2) {
		jags.rem 		<- which(jags.params == "rho.mat")
		jags.params	<- jags.params[-jags.rem]
	}
	#function used to initialize parameters
    jags.inits <- list( dmu.prior.mu=dmu.prior.mu, mu.prior.mu=mu.prior.mu, p.transform=runif(num.sources), region.sig=0.5, ind.sig=0.5, p.ind = matrix(runif(N*num.sources), N, num.sources) )
# print(curr.model)
# print(file.flag)
	if(run.parallel) {parallel.state <- "parallel"} else { parallel.state <- "interruptible"}
	jags.out <- run.jags(model=curr.model, monitor=jags.params, data=jags.dump, n.chains=mcmc.chains, burnin=mcmc.burn, sample= (mcmc.chainLength-mcmc.burn), thin=mcmc.thin, psrf.target=1.5, check.conv=TRUE, plots=FALSE, check.stochastic=FALSE,  monitor.deviance=TRUE, silent.jags=FALSE, monitor.pd=FALSE, method=parallel.state)

	r.est <- jags.out$psrf$psrf[,1]
	jags.output.mat <- cbind(jags.out$summary$statistics[,1:2], jags.out$summary$quantiles, r.est)
    
	save(jags.out, jags.output.mat, X, sources, nome.flag, num.sources, num.iso, mcmc.chains, N, num.groups, file=output.name)
    write.table(jags.output.mat, file=paste(strsplit(output.name, ".Rdata")[[1]],'.txt',sep=''))
	

    if(plot.observations) {
		if(num.sources >= 3 & num.iso==2) {dev.new(); Tri.plots(jags.out, X, sources=sources, plot.ind=TRUE, me.flag=!nome.flag, color.plots=color.plots) } else {
			if(num.sources == 2 & num.iso==2) {dev.new(); Bi.plots(jags.out, X, sources=sources, plot.ind=TRUE, me.flag=!nome.flag, color.plots=color.plots)} else {warning("No observation plot available for this number of isotopes", call.=FALSE) }
		}
    }
    options(warn=0)

    if(plot.mixing.estimates) {
        
        if(num.sources == 3 & num.iso==2) {dev.new(); Tri.plots(jags.1=jags.out, X=X, sources=sources, plot.mix=TRUE, me.flag=!nome.flag, color.plots=color.plots) } else {
        if(num.sources == 2 & num.iso==2) {dev.new(); Bi.plots(jags.out, X,  sources=sources, plot.mix=TRUE, me.flag=!nome.flag, color.plots=color.plots) } else {warning("No mixing space plot available for isotope/source combination", call.=FALSE)}
        }    
    }
    
    if(plot.dietary.source.contributions) {
        dev.new()
        curves.plot(jags.1=jags.out, num.sources=num.sources, num.chains=mcmc.chains, color=color.plots, individuals=N,  xlab.vec=levels(as.factor(sources[,num.iso+1])), num.groups=num.groups)
    }
	print(jags.output.mat)

}

##load previously run data and outputs graphs
load.prev.func <- function(file.name="SampleOutput.Rdata", plot.observations=TRUE, plot.mixing.estimates=TRUE, plot.dietary.source.contributions=TRUE, color.plots=TRUE) {

	jags.out=NA
# 	X=NA
	sources=NA
	nome.flag=NA
# 	num.iso=NA
	num.sources=NA
	mcmc.chains=NA
	N=NA
	num.groups=NA
	jags.output.mat=NA
	load(file=file.name)

	if(!exists("jags.out") | !exists("X") | !exists("sources") | !exists("nome.flag") | !exists("num.iso") | !exists("num.sources") | !exists("mcmc.chains") | !exists("N") | !exists("num.groups") | !exists("jags.output.mat"))
	{stop(".Rdata file error")}
	
	if(plot.observations) {
		if(num.sources >= 3 & num.iso==2) {dev.new(); Tri.plots(jags.out, X, sources=sources, plot.ind=TRUE, me.flag=!nome.flag, color.plots=color.plots) } else {
			if(num.sources == 2 & num.iso==2) {dev.new(); Bi.plots(jags.out, X, sources=sources, plot.ind=TRUE, me.flag=!nome.flag, color.plots=color.plots)} else {warning("No observation plot available for this number of isotopes", call.=FALSE) }
		}

		# 		  if(num.iso==2) {dev.new(); Tri.plots(jags.out, X, sources=sources, plot.ind=TRUE, me.flag=!nome.flag, color.plots=color.plots) } 
    }

    if(plot.mixing.estimates) {
        if(num.iso == 2 & num.sources==3) {dev.new(); Tri.plots(jags.1=jags.out, X=X, sources=sources, plot.mix=TRUE, me.flag=!nome.flag) } else {
        if(num.iso == 2 & num.sources==2) {dev.new(); Bi.plots(jags.out, X,  sources=sources, plot.mix=TRUE, me.flag=!nome.flag) } else {warning("Warning: No mixing plot available for this isotope/source combination", call.=FALSE)}
        }    
    }
    
 
    if(plot.dietary.source.contributions) {
        dev.new()
        curves.plot(jags.1=jags.out, num.sources=num.sources, num.chains=mcmc.chains, color=color.plots, individuals=N,  xlab.vec=levels(as.factor(sources[,num.iso+1])), num.groups=num.groups)
    }
    print(jags.output.mat)
	return(NULL)
}

IsotopeR    <- function() {
    require(fgui) #used to run JAGS from R

    if(interactive()) {
		fguiWindow( basicMenu=FALSE, title="IsotopeR", text="Please choose an option from the Analysis menu." )
		win <- mgui(IsoWrapper, argFilter = list(Mixtures="{{} {.csv}}", Sources="{{} {.csv}}", Concentrations="{{} {.csv}}", Measurement.Error="{{} {.csv}}",   Discrimination.Error="{{} {.csv}}", Digestibility.Factor="{{} {.csv}}"), argText = list(mcmc.chains="number of chains", mcmc.burn="MCMC burnin", mcmc.chainLength="MCMC runs", mcmc.thin="thinning rate", output.name="Output file"), argOption = list(run.parallel=c("TRUE", "FALSE"), plot.observations=c("TRUE","FALSE"), plot.dietary.source.contributions= c("TRUE", "FALSE"), plot.mixing.estimates=c('TRUE','FALSE'), color.plots=c("TRUE", "FALSE")), closeOnExec=TRUE, title=c("Analysis","New Run"), exec="Run IsotopeR", helpsFunc="IsoWrapper", output=NULL)
		
		win <- mgui(load.prev.func, title=c("Analysis","Load Previous Run"), argFilter= list(file.name="{{} {.Rdata}}"), argOption = list(plot.observations=c("TRUE","FALSE"), plot.dietary.source.contributions= c("TRUE", "FALSE"), plot.mixing.estimates=c('TRUE','FALSE'), color.plots=c("TRUE", "FALSE")), closeOnExec=TRUE , exec="Plot", output=NULL, helpsFunc="IsoWrapper")
# 		print(win)
		
    }


}

