##Isotope Mixing ,Model with Measurement Error, Source Correlation, Concentration Error and residual error
##doesnt estimate third source concentration, inputs mean and variance
##outputs C:N ratios
##Jake Ferguson updated 5/1/2011
##based on original code from moore and semmens 2009

IsotopeRfull <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
  ##covariance matrix
  for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    
    for(source2 in 1:num.iso) {
      D.tau.temp[source2,source2,source] ~ dgamma(.001,.001)
    }
  
    D.tau[1:num.iso,1:num.iso,source] <- D.tau.temp[,,source]
#     mu.conc[1:num.iso,source] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau[,])
      
#     for(count in 1:cd.ss[source]) {
#         cd.array[count,1:num.iso,source] ~ dmnorm(mu.conc[,source], D.tau[,,source])
#     }


	##draws subsource means and fits to data
    for(sub in 1:subcd.vec[source]) { 
		subcd[source, sub, 1:num.iso] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau)
	}
	
	for(iso in 1:num.iso) {
		mu.conc[iso, source] <- mean( subcd[source, 1:subcd.vec[source], iso])
	}
		
	for(sub in 1:subcd.vec[source]) {
		for(index  in subcd.samples[source,sub,1]:subcd.samples[source,sub,2]) {
			cd.mat[index, ] ~ dmnorm(subcd[source,sub, ], D.tau[ , ,source]);
		}  
     }

  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.1,.1)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.1,.1)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  cov.ME + discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }

  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source]+discrimvar.mat[iso,iso])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {
        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
        sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])
    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconc <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  cov.ME + discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }

	}#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source]+discrimvar.mat[iso,iso])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconcnodiscrim <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
   for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  cov.ME )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconc <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  cov.ME + discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source]+discrimvar.mat[iso,iso])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconcnodiscrim <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  cov.ME )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconcnome <- "model {
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
      for(iso in 1:num.iso) {
        mu.conc[iso,source] <- 1
      }
  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source])  + discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source]+discrimvar.mat[iso,iso])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
#     sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
#         sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnoconcnomenodiscrim <-  "model {

	for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        mu.conc[iso,source] <- 1
      }
  }
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.01);} 
  pop.invSig2 ~ dgamma(.01, .01);
  pop.var <- 1/pop.invSig2;
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source], pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01, .01);
  ind.var <- 1/ind.invSig2;
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop;
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso];
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
		p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,];
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso];
      }
    }
  }#end for i
  
  ##################
  ##Source Estimation####
  ##################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        tau.source.temp[sourcex,sourcey,source] <- 0;
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        tau.source.temp[sourcex,sourcey,source] <- 0;
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.001, .001);
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1);
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source];
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1;
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source]);
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source])  ); #+ discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
      
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(.01, .01);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau);

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source];
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0;
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0;
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  res.err);
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
        sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i]);
      }
    }
    
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso]);
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso]);
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source]);

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind] );
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }
}"

IsotopeRnodiscrim   <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
  for(iso1 in 2:num.iso) {
    for(iso2 in 1:(iso1-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso2 in 2:num.iso) {
    for(iso1 in 1:(iso2-1)) {
        tauZ[iso1,iso2] <- 0
    }
  }
  for(iso in 1:num.iso) {
    tauZ[iso,iso] ~ dgamma(.001,.001)
  }

  cov.ME <- inverse(tauZ);

  for(i in 1:Nz) {
    Z[i,] ~ dmnorm(muz,tauZ);
  }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
  ##covariance matrix
  for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    
    for(source2 in 1:num.iso) {
      D.tau.temp[source2,source2,source] ~ dgamma(.001,.001)
    }
  
    D.tau[1:num.iso,1:num.iso,source] <- D.tau.temp[,,source]
#     mu.conc[1:num.iso,source] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau[,])
      
#     for(count in 1:cd.ss[source]) {
#         cd.array[count,1:num.iso,source] ~ dmnorm(mu.conc[,source], D.tau[,,source])
#     }
      
	##draws subsource means and fits to data
    for(sub in 1:subcd.vec[source]) { 
		subcd[source, sub, 1:num.iso] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau)
	}
	
	for(iso in 1:num.iso) {
		mu.conc[iso, source] <- mean( subcd[source, 1:subcd.vec[source], iso])
	}
		
	for(sub in 1:subcd.vec[source]) {
		for(index  in subcd.samples[source,sub,1]:subcd.samples[source,sub,2]) {
			cd.mat[index, ] ~ dmnorm(subcd[source,sub, ], D.tau[ , ,source]);
		}  
     }

  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.1,.1)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.1,.1)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  cov.ME )#+ discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
    for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
	for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  cov.ME + res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
    sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {
        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
        sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])
    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


IsotopeRnome <- "model {

  ##################################
  ##estimate the measurement error##
  ##################################
#   for(iso in 1:num.iso) {
#     tauZ[iso,iso] ~ dgamma(.001,.001)
#   }
 
#   cov.ME <- 0#inverse(tauZ);

#   for(i in 1:Nz) {
#     Z[i,] ~ dmnorm(muz,tauZ);
#   }
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
  ##covariance matrix
  for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    
    for(source2 in 1:num.iso) {
      D.tau.temp[source2,source2,source] ~ dgamma(.001,.001)
    }
  
    D.tau[1:num.iso,1:num.iso,source] <- D.tau.temp[,,source]
#     mu.conc[1:num.iso,source] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau[,])
      
#     for(count in 1:cd.ss[source]) {
#         cd.array[count,1:num.iso,source] ~ dmnorm(mu.conc[,source], D.tau[,,source])
#     }
      
	##draws subsource means and fits to data
    for(sub in 1:subcd.vec[source]) { 
		subcd[source, sub, 1:num.iso] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau)
	}
	
	for(iso in 1:num.iso) {
		mu.conc[iso, source] <- mean( subcd[source, 1:subcd.vec[source], iso])
	}
		
	for(sub in 1:subcd.vec[source]) {
		for(index  in subcd.samples[source,sub,1]:subcd.samples[source,sub,2]) {
			cd.mat[index, ] ~ dmnorm(subcd[source,sub, ], D.tau[ , ,source]);
		}  
     }

  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
    p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
    tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
      tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]) +  discrimvar.mat )
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source]+discrimvar.mat[iso,iso])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
#     sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
        sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}"


##debuggin this one
IsotopeRnomenodiscrim <- "model {
    
  ###############################
  ##estimate the concentrations##
  ###############################
  for(source in 1:(num.sources)) {
  ##covariance matrix
  for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
        D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        D.tau.temp[sourcex,sourcey,source] <- 0
      }
    }    
    for(source2 in 1:num.iso) {
      D.tau.temp[source2,source2,source] ~ dgamma(.001,.001)
    }
  
    D.tau[1:num.iso,1:num.iso,source] <- D.tau.temp[,,source]
#     mu.conc[1:num.iso,source] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau[,])
      
#     for(count in 1:cd.ss[source]) {
#         cd.array[count,1:num.iso,source] ~ dmnorm(mu.conc[,source], D.tau[,,source])
#     }      


	##draws subsource means and fits to data
    for(sub in 1:subcd.vec[source]) { 
		subcd[source, sub, 1:num.iso] ~ dmnorm(dmu.prior.mu[,source], dmu.prior.tau)
	}
	
	for(iso in 1:num.iso) {
		mu.conc[iso, source] <- mean( subcd[source, 1:subcd.vec[source], iso])
	}
		
	for(sub in 1:subcd.vec[source]) {
		for(index  in subcd.samples[source,sub,1]:subcd.samples[source,sub,2]) {
			cd.mat[index, ] ~ dmnorm(subcd[source,sub, ], D.tau[ , ,source]);
		}  
     }

  }
  
  
  ###########################
  ##Proportion estimamation##
  ###########################
  ## this is the global mean
  for(i in 1:num.sources) {mu[i] ~ dnorm(alpha.clr[i], 0.001)} 
  pop.invSig2 ~ dgamma(.01,.01)
  pop.var <- 1/pop.invSig2
  for(source in 1:num.sources) {
    p.transform[source] ~ dnorm(mu[source],pop.invSig2);
  }  
   
  ##generate individuals draws from the global mean
  ind.invSig2 ~ dgamma(.01,.01)
  ind.var <- 1/ind.invSig2
  for(i in 1:N) {
    for(source in 1:num.sources) {
      p.ind[i,source] ~ dnorm(p.transform[source], ind.invSig2);
      exp.p[i,source] <- exp(p.ind[i,source]);
    }
  }
    
  ##CLR math: This does the back-transform to get back to proportions
  for(source in 1:num.sources) {
    p.exp.transform[source] <- exp(p.transform[source]);
  }
  p.mean.tot <- sum(p.exp.transform[]);
  for(source in 1:(num.sources-1)) {
    p.pop[source] <- p.exp.transform[source]/p.mean.tot;
  }
  p.pop[num.sources] <- 1-sum(p.pop[1:(num.sources-1)]); 
  
  ##rescale p.pop for concentration dependence
  p.popdenom <- mu.conc%*%p.pop
  for(iso in 1:num.iso) {
    for(source in 1:num.sources) {
      pIso.pop[iso,source]  <- mu.conc[iso,source]*p.pop[source]/p.popdenom[iso]
    }
  }
  
  ##individual p's
  for(i in 1:N) {
    p.tot[i] <- sum(exp.p[i,1:num.sources]);
      for(source in 1:(num.sources-1)) {
        p[i,source] <- exp.p[i,source]/p.tot[i];
      }
    p[i,num.sources] <- 1-sum(p[i,1:(num.sources-1)]);
   
    ##rescale p.pop for concentration dependence
    p.inddenom[i,1:num.iso] <- mu.conc%*%p[i,]
    for(iso in 1:num.iso) {
      for(source in 1:num.sources) {
        pIso.ind[i,iso,source]  <- mu.conc[iso,source]*p[i,source]/p.inddenom[i,iso]
      }
    }
 }#end for i
  
  #####################
  ##Source Estimation##
  #####################
  ##estimate sources
  for(source in 1:num.sources) { 
   ##covariance matrix
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
			tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
		tau.source.temp[sourcex,sourcey,source] <- 0
      }
    }
    for(source2 in 1:num.iso) {
		tau.source.temp[source2,source2,source] ~ dgamma(.0010,.0010)
    }
       
    ##build source correlation matrix
    rho.source[source] ~ dunif(-1,1)
    for(sourcex in 2:(num.iso)) {
      for(sourcey in 1:(sourcex-1)) {
		rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(sourcey in 2:(num.iso)) {
      for(sourcex in 1:(sourcey-1)) {
        rho.mat[sourcex,sourcey,source] <- rho.source[source]
      }
    }
    for(source2 in 1:num.iso) {
      rho.mat[source2,source2,source] <- 1 
    }      
    cov.source[1:num.iso,1:num.iso,source] <- inverse(tau.source.temp[,,source])
    tau.source[1:num.iso,1:num.iso,source] <- inverse((cov.source[,,source]%*%rho.mat[,,source]%*%cov.source[,,source]))
 
    ##draws subsource means and fit to data
    for(sub in 1:subsource.vec[source]) { 
		subsource[source, sub, 1:num.iso] ~ dmnorm(mu.prior.mu, mu.prior.cov)
	}
	
	for(iso in 1:num.iso) {
		mu.source[source, iso] <- mean( subsource[source, 1:subsource.vec[source], iso])
	}
		
	for(sub in 1:subsource.vec[source]) {
		for(index  in subsource.samples[source,sub,1]:subsource.samples[source,sub,2]) {
			source.mat[index, ] ~ dmnorm(subsource[source,sub, ], tau.source[ , ,source]);
		}  
     }
  }#end for sources
  
  
  ####################### 
  ##draw residual error##
  #######################
  for(i in 1:(num.iso-1)) {
	  for(j in (i+1):num.iso) {
		res.tau[i,j] <- 0;
		res.tau[j,i] <- 0;
		}
	}
  for(iso in 1:num.iso) {
    res.tau[iso,iso] ~ dgamma(1e-3,1e-3)#dexp(1/1000)#dunif(0,20)#dgamma(10,10)#dunif(0,20);#dexp(1);
  }
  res.err[1:num.iso,1:num.iso] <- inverse(res.tau)

  ##rescale sources by p
  for(i in 1:N) {
    ##rescale covariance and include fractionation
    for(source in 1:num.sources) {
      for(iso in 1:num.iso) {
        covfrac.source[iso,iso,source,i] <- (cov.source[iso,iso,source])*pIso.ind[i,iso,source]
      }
      for(isox in 2:(num.iso)) {
    for(isoy in 1:(isox-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
      for(isoy in 2:(num.iso)) {
    for(isox in 1:(isoy-1)) {
      covfrac.source[isox,isoy,source,i] <- 0
    }
      }
  
      obscov.mat[1:num.iso,1:num.iso,source,i] <- (covfrac.source[,,source,i]%*%rho.mat[,,source]%*%covfrac.source[,,source,i] +  res.err)  
  
    }#end for source

    for(x in 1:num.iso) {
      for(y in 1:num.iso) {
    sumobscov.mat[x,y,i] <- sum(obscov.mat[x,y,1:num.sources,i])  
      }
    }
    for(iso in 1:num.iso) {
      mu.mix[i,iso] <- pIso.ind[i,iso,]%*%(mu.source[,iso])
    }
  }#end for i

  ##get the sd's for jack
  for(iso in 1:num.iso) {
  
    sd.res[iso] <- sqrt(res.err[iso,iso])
#     sd.me[iso] <- sqrt(cov.ME[iso,iso])
            
    for(source in 1:num.sources) {

        sd.source[iso,source] <- sqrt(cov.source[iso,iso,source])
        sd.conc[iso,source] <- 1/sqrt(D.tau[iso,iso,source])

    }
  }


  ##calculate the likelihoods for the N individuals.
  for(ind in 1:N) {
    mix.prcsn[1:num.iso,1:num.iso,ind] <- inverse(sumobscov.mat[,,ind]  )
    for(j in 1:ind.counts[ind]) {      
      ind.array[1:num.iso,ind,j] ~ dmnorm(mu.mix[ind,1:num.iso], mix.prcsn[1:num.iso,1:num.iso, ind]);    
    }
  }

}#end model
"

