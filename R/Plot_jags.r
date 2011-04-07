########output plots for methods paper###########
####Jake Ferguson 1/7/11

##################################
##Triangle plot in isotope space
##jags.1 is an jags object from an IsotopeR model run
##X is the observed isotope mixture values for each individual
##plot.mix is a 0,1 flag denoting whether to plot estimated mixture values
##plot.ind.flag is a 0,1 flag denoting whether to plot X.
##################################
Tri.plots <- function(jags.1,X,plot.mix=FALSE,plot.ind.flag=FALSE, me.flag=FALSE) {
  require(ellipse)
  require(plotrix)
  
  num.ind <- dim(X)[1]

  source.cols <- c("black","black","black")#c("chartreuse3", "goldenrod2","cornflowerblue")#,)#multiline.plot.colors()#c
  jags.names <- dimnames(jags.1$summary$quantiles)
  jags.table <- jags.1$summary$quantiles
  med <- which(jags.names[[2]] == '50%')
  low.ci <- which(jags.names[[2]] == '2.5%')
  hi.ci <- which(jags.names[[2]] == '97.5%')
  
  mu.source.index <- grep("mu.source", jags.names[[1]])
  mu.source <- jags.table[mu.source.index,med]
  mu.lowCI <- jags.table[mu.source.index,low.ci]
  mu.hiCI <- jags.table[mu.source.index,hi.ci]
  cov.source.index <- grep("sd.source", jags.names[[1]])
  cov.source <- jags.table[cov.source.index, med]

  if(me.flag) {
    sigmaz.index <- grep("sd.me", jags.names[[1]])
    sigmaz.temp <- jags.table[sigmaz.index, med]
    sigmaz.med  <- sigmaz.temp#[c(1,4)]
  }
  
  mix.index <- grep("mu.mix", jags.names[[1]])
  mix.temp <- jags.table[mix.index,med]
  mixmu.med 	<- cbind(mix.temp[1:num.ind], mix.temp[(num.ind+1):(2*num.ind)])

  mixmu.loCI 	<- jags.table[mix.index,low.ci]
  mixmu.loCI 	<- cbind(mixmu.loCI[1:num.ind], mixmu.loCI[(num.ind+1):(2*num.ind)])
  mixmu.hiCI	<- jags.table[mix.index,hi.ci]
  mixmu.hiCI 	<- cbind(mixmu.hiCI[1:num.ind], mixmu.hiCI[(num.ind+1):(2*num.ind)])
  
  cov.conc.index <- grep("sd.conc", jags.names[[1]])
  cov.conc <- jags.table[cov.conc.index, med]
  
  x.points <- mu.source[1:(length(mu.source)/2)]
  if(me.flag) { x.sd     <- cov.source[1:3] + sigmaz.med[1] } else { 
    x.sd		<- cov.source[1:3] #+ sigmaz.med[1]^2)
  }
  x.lowCI	<- x.points - 2*x.sd
  x.hiCI	<- x.points + 2*x.sd

  y.points <- mu.source[(length(mu.source)/2+1):length(mu.source)]
  y.lowCI <- mu.lowCI[(length(mu.source)/2+1):length(mu.source)]
  y.hiCI <- mu.hiCI[(length(mu.source)/2+1):length(mu.source)]
	
  if(me.flag) {
      y.sd <- cov.source[4:6]+ sigmaz.med[2]
  } else {
      y.sd <- cov.source[4:6]
  }
  y.lowCI <- y.points - 2*y.sd
  y.hiCI <- y.points + 2*y.sd

  data.cols="black"
  
  D.source.index <- grep("mu.conc", jags.names[[1]])
  D.source <- jags.table[D.source.index,med]

  D.lowCI <- jags.table[D.source.index,low.ci]
  D.hiCI <- jags.table[D.source.index,hi.ci]
  
  x.coord <- seq(1,length(D.source),by=2)
  x.conc <- D.source[x.coord]
  x.conc.sd <- cov.source[1:3] #sqrt(cov.source[seq(1,12,by=4)])
 
  x.conc.lowCI <- x.conc-2*x.conc.sd
  x.conc.hiCI <- x.conc-2*x.conc.sd

  y.coord <- seq(2,length(D.source),by=2)
  y.conc <- D.source[y.coord]
  y.conc.sd <- cov.source[4:6]

  y.conc.lowCI <- y.conc - 2*y.conc.sd
  y.conc.hiCI <- y.conc - 2*y.conc.sd
  
  axis.raw <- seq(0,1,by=0.2)
  base.matrix <- array(NA,c(length(axis.raw),length(axis.raw),3))
  x.matrix <- base.matrix
  xedge.matrix <- base.matrix
  x.lowCI.matrix <- base.matrix
  x.hiCI.matrix <- base.matrix
  y.matrix <- base.matrix
  yedge.matrix <- base.matrix
  y.lowCI.matrix <- base.matrix
  y.hiCI.matrix <- base.matrix
  x.outer.matrix <- base.matrix
  y.outer.matrix <- base.matrix
  
  for(r1 in 1:length(axis.raw)) {
    for(r2 in 1:(length(axis.raw))) {
      temp <- 1- axis.raw[r1] - axis.raw[r2] 
      if(temp >= (0.00-1e-3) & temp <=1.001) {
	base.matrix[r1,r2,1] <- axis.raw[r1]
	base.matrix[r1,r2,2] <- axis.raw[r2]
	base.matrix[r1,r2,3] <- 1- axis.raw[r1] - axis.raw[r2]
	
	x.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*x.conc[1]/base.matrix[r1,r2,]%*%x.conc
	x.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*x.conc[2]/base.matrix[r1,r2,]%*%x.conc
	x.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*x.conc[3]/base.matrix[r1,r2,]%*%x.conc
	
	y.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*y.conc[1]/base.matrix[r1,r2,]%*%y.conc
	y.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*y.conc[2]/base.matrix[r1,r2,]%*%y.conc
	y.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*y.conc[3]/base.matrix[r1,r2,]%*%y.conc

	x.lowCI.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*x.conc.lowCI[1]/base.matrix[r1,r2,]%*%x.conc.lowCI
	x.lowCI.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*x.conc.lowCI[2]/base.matrix[r1,r2,]%*%x.conc.lowCI
	x.lowCI.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*x.conc.lowCI[3]/base.matrix[r1,r2,]%*%x.conc.lowCI

	y.lowCI.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*y.conc.lowCI[1]/base.matrix[r1,r2,]%*%y.conc.lowCI
	y.lowCI.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*y.conc.lowCI[2]/base.matrix[r1,r2,]%*%y.conc.lowCI
	y.lowCI.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*y.conc.lowCI[3]/base.matrix[r1,r2,]%*%y.conc.lowCI

	x.hiCI.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*x.conc.hiCI[1]/base.matrix[r1,r2,]%*%x.conc.hiCI
	x.hiCI.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*x.conc.hiCI[2]/base.matrix[r1,r2,]%*%x.conc.hiCI
	x.hiCI.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*x.conc.hiCI[3]/base.matrix[r1,r2,]%*%x.conc.hiCI

	y.hiCI.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*y.conc.hiCI[1]/base.matrix[r1,r2,]%*%y.conc.hiCI
	y.hiCI.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*y.conc.hiCI[2]/base.matrix[r1,r2,]%*%y.conc.hiCI
	y.hiCI.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*y.conc.hiCI[3]/base.matrix[r1,r2,]%*%y.conc.hiCI
	

	if(r1==1 | r2==1 | r1+r2 == (length(axis.raw)+1)) {
	  xedge.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*x.conc[1]/base.matrix[r1,r2,]%*%x.conc
	  xedge.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*x.conc[2]/base.matrix[r1,r2,]%*%x.conc
	  xedge.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*x.conc[3]/base.matrix[r1,r2,]%*%x.conc
	  
	  yedge.matrix[r1,r2,1] <- base.matrix[r1,r2,1]*y.conc[1]/base.matrix[r1,r2,]%*%y.conc
	  yedge.matrix[r1,r2,2] <- base.matrix[r1,r2,2]*y.conc[2]/base.matrix[r1,r2,]%*%y.conc
	  yedge.matrix[r1,r2,3] <- base.matrix[r1,r2,3]*y.conc[3]/base.matrix[r1,r2,]%*%y.conc
	  
	}	
      }
    }
  }#end r loop

  x.basic <- apply(x.matrix, c(1,2), '%*%', x.points)
  y.basic <- apply(y.matrix, c(1,2), '%*%', y.points)
  
  x.edge <- apply(xedge.matrix, c(1,2), '%*%', x.points)
  y.edge <- apply(yedge.matrix, c(1,2), '%*%', y.points)

  x.edge.lowCI <- apply(x.lowCI.matrix, c(1,2), '%*%', x.points)
  y.edge.lowCI <- apply(y.lowCI.matrix, c(1,2), '%*%', y.points)

  x.edge.hiCI <- apply(x.hiCI.matrix, c(1,2), '%*%', x.points)
  y.edge.hiCI <- apply(y.hiCI.matrix, c(1,2), '%*%', y.points)
  
  x.outer.loCI <- apply(x.matrix, c(1,2), '%*%', x.lowCI)
  x.outer.hiCI <- apply(x.matrix, c(1,2), '%*%', x.hiCI)

  y.outer.loCI <- apply(y.matrix, c(1,2), '%*%', y.lowCI)
  y.outer.hiCI <- apply(y.matrix, c(1,2), '%*%', y.hiCI)

  plotCI(x=x.points, y=y.points, liw=(x.sd), uiw=(x.sd), xlim=(range(x.points) + c(-1,1)*2*max(x.sd)), ylim=(range(y.points) + c(-1,1)*2*max(y.sd)), err="x",pch=19, xlab=expression(paste(delta^13,C)) , ylab=expression(paste(delta^15,N)),col="white")
  plotCI(x=x.points, y=y.points, liw=(y.sd), uiw=(y.sd) ,err="y",add=TRUE,pch=19,col="white")
  box(lwd=2)
    
  base.matrix <- round(base.matrix,1)
  ##draws the grey triangles
  for(x1 in 1:(length(axis.raw))) {
    for(y1 in 1:(length(axis.raw))) {
      
    for(x2 in (x1-1):(x1+1)) {
	for(y2 in (y1-1):(y1+1)) {
	  
	  if(x2 > (length(axis.raw)) | x2 < 1.0 | y2 > (length(axis.raw)) | y2 < 1.0 ) {next()} else { 

	    segments(x0=x.basic[x1,y1], y0=y.basic[x1,y1], x1=x.basic[x2,y2], y1=y.basic[x2,y2],col="grey",lwd=1,lty=1)
	  }
	}
      }
    }#end for x1
  }#end for y1
  
  ##draws the black edges & the CIs
  edge.points <- which(!is.na(x.edge))
  edge2 <- c(edge.points,edge.points)
  back.vec <- (length(axis.raw)):1

  for(i in 1:(length(axis.raw)-1)) {

    ##middle
    segments(x0=x.basic[i,1], y0=y.basic[i,1], x1=x.basic[i+1,1], y1=y.basic[i+1,1],col="darkgrey",lwd=1.5)
    segments(x0=x.basic[1,i], y0=y.basic[1,i], x1=x.basic[1,i+1], y1=y.basic[1,i+1],col="darkgrey",lwd=1.5)
    segments(x0=x.basic[back.vec[i],i], y0=y.basic[back.vec[i],i], x1=x.basic[back.vec[i]-1,i+1], y1=y.basic[back.vec[i]-1,i+1],col="darkgrey",lwd=1.5)
    
  }
  
  points(x.points,y.points,pch=19,col="grey")
  points(x.points,y.points,pch=15,col=source.cols,lwd=2)
  theta <- seq(0, 2 * pi, length=100)

  rho.index <- grep("rho.source", jags.names[[1]])
  rho.vec <-jags.1$summary$statistics[rho.index,1]
 
  for(i in 1:3) {

    T.mat <- diag(2)
    T.mat[1,2] <- rho.vec[i]
    T.mat[2,1] <- rho.vec[i]
    
    sd.vec <- c(cov.source[i], cov.source[i+3])
    if(me.flag) {me.mat <- diag(sigmaz.med[1:2]^2)
         T.mat  <- c(cov.source[i], cov.source[i], cov.source[i+3], cov.source[i+3])*T.mat*sd.vec + me.mat
    } else {
         T.mat  <- c(cov.source[i], cov.source[i], cov.source[i+3], cov.source[i+3])*T.mat*sd.vec
    }

    lines(ellipse(T.mat, centre=c(x.points[i], y.points[i]), level=0.95), lty=2, lwd=2)
  }
  
  
  ##plots predicted isotope values of individuals
  if(plot.mix ==TRUE) {
    plotCI(x=mixmu.med[,1], y=mixmu.med[,2], liw=(mixmu.med[,1]-mixmu.loCI[,1]), uiw=(mixmu.hiCI[,1] - mixmu.med[,1]), err="x", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col='dimgrey', pch=19)
    plotCI(x=mixmu.med[,1], y=mixmu.med[,2], liw=(mixmu.med[,2]-mixmu.loCI[,2]), uiw=(mixmu.hiCI[,2] - mixmu.med[,2]), err="y", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col='dimgrey', pch=19)
    points(mixmu.med,lwd=2)
    title("Estimates")
  }

 ##Plots observed istope values of individuals
    if(plot.ind.flag) {  
        if(me.flag) {
            plotCI(x=X[,1], y=X[,2], uiw = 2*sqrt(sigmaz.med[1]), err="x", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col='dimgrey', pch=19)
            plotCI(x=X[,1], y=X[,2], uiw = 2*sqrt(sigmaz.med[2]), err="y", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col='dimgrey', pch=19)
        }
        points(X,lwd=2)
        title("Observations")
    }  
  
}#end triplots



##################################
##2 source plot in isotope space
##jags.1 is an jags object from an IsotopeR model run
##X is the observed isotope mixture values for each individual
##plot.mix is a 0,1 flag denoting whether to plot estimated mixture values
##plot.ind.flag is a 0,1 flag denoting whether to plot X
##################################
Bi.plots <- function(jags.1,X,plot.mix=FALSE,plot.ind.flag=FALSE, me.flag=FALSE) {
  
  require(ellipse)
  require(plotrix)
  
  N <- dim(X)[1]
  
  jags.names <- dimnames(jags.1$BUGSoutput$summary)
  jags.table <- jags.1$BUGSoutput$summary
  med <- which(jags.names[[2]] == '50%')
  low.ci <- which(jags.names[[2]] == '2.5%')
  hi.ci <- which(jags.names[[2]] == '97.5%')
  
  mu.source.index <- grep("mu.source", jags.names[[1]])
  mu.source <- jags.table[mu.source.index,med]
  mu.lowCI <- jags.table[mu.source.index,low.ci]
  mu.hiCI <- jags.table[mu.source.index,hi.ci]
 
  cov.source.index <- grep("sd.source", jags.names[[1]])
  cov.source <- jags.table[cov.source.index, med]

    if(me.flag) {
        sigmaz.index <- grep("sd.me", jags.names[[1]])
        sigmaz.temp <- jags.table[sigmaz.index, med]
        sigmaz.med  <- sigmaz.temp#[c(1,4)]
    }

  mix.index <- grep("mu.mix", jags.names[[1]])
  mix.temp <- jags.table[mix.index,med]
  
  mixmu.med 	<- cbind(mix.temp[1:N], mix.temp[(N+1):(2*N)])
  mixmu.loCI 	<- jags.table[mix.index,low.ci]
  mixmu.loCI 	<- cbind(mixmu.loCI[1:N], mixmu.loCI[(N+1):(2*N)])
  mixmu.hiCI	<- jags.table[mix.index,hi.ci]
  mixmu.hiCI 	<- cbind(mixmu.hiCI[1:N], mixmu.hiCI[(N+1):(2*N)])
  
  cov.conc.index <- grep("sd.conc", jags.names[[1]])
  cov.conc <- jags.table[cov.conc.index, med]
  
  x.points <- mu.source[1:(length(mu.source)/2)]

  if(me.flag) { x.sd     <- cov.source[1:2] + sigmaz.med[1] } else { 
    x.sd        <- cov.source[1:2] #+ sigmaz.med[1]^2)
  }
  x.lowCI <- x.points-2*x.sd
  x.hiCI <- x.points+2*x.sd
  
  y.points <- mu.source[(length(mu.source)/2+1):length(mu.source)]
  y.lowCI <- mu.lowCI[(length(mu.source)/2+1):length(mu.source)]
  y.hiCI <- mu.hiCI[(length(mu.source)/2+1):length(mu.source)]

  if(me.flag) {
      y.sd <- cov.source[3:4]+ sigmaz.med[2]
  } else {
      y.sd <- cov.source[3:4]
  }
  y.lowCI <- y.points - 2*y.sd
  y.hiCI <- y.points + 2*y.sd
  
  data.cols="black"
  
  D.source.index <- grep("mu.conc", jags.names[[1]])
  D.source <- jags.table[D.source.index,med]

  D.lowCI <- jags.table[D.source.index,low.ci]
  D.hiCI <- jags.table[D.source.index,hi.ci]
  
  x.coord <- seq(1,length(D.source),by=2)
  x.conc <- D.source[x.coord]
	x.conc.sd <- cov.source[1:2]#sqrt(cov.source[seq(1,8,by=4)])

  x.conc.lowCI <- x.conc-2*x.conc.sd
  x.conc.hiCI <- x.conc+2*x.conc.sd
  
  y.coord <- seq(2,length(D.source),by=2)
  y.conc <- D.source[y.coord]
  y.conc.sd <- cov.source[3:4]#sqrt(cov.source[seq(4,8,by=4)])
  
  y.conc.lowCI <- y.conc-2*y.conc.sd
  y.conc.hiCI <- y.conc+2*y.conc.sd
  
  axis.raw <- seq(0,1,by=0.2)
  base.matrix <- array(NA,c(length(axis.raw),2))
  x.matrix <- base.matrix
  xedge.matrix <- base.matrix
  x.lowCI.matrix <- base.matrix
  x.hiCI.matrix <- base.matrix
  y.matrix <- base.matrix
  yedge.matrix <- base.matrix
  y.lowCI.matrix <- base.matrix
  y.hiCI.matrix <- base.matrix
  x.outer.matrix <- base.matrix
  y.outer.matrix <- base.matrix
  
  for(r1 in 1:length(axis.raw)) {
      temp <- 1- axis.raw[r1] #- axis.raw[r2] 
      if(temp >= (0.00-1e-3) & temp <=1.001) {
	base.matrix[r1,1] <- axis.raw[r1]
	base.matrix[r1,2] <- 1-axis.raw[r1]#axis.raw[r2]
	
	x.matrix[r1,1] <- base.matrix[r1,1]*x.conc[1]/base.matrix[r1,]%*%x.conc
	x.matrix[r1,2] <- base.matrix[r1,2]*x.conc[2]/base.matrix[r1,]%*%x.conc

	y.matrix[r1,1] <- base.matrix[r1,1]*y.conc[1]/base.matrix[r1,]%*%y.conc
	y.matrix[r1,2] <- base.matrix[r1,2]*y.conc[2]/base.matrix[r1,]%*%y.conc

	x.lowCI.matrix[r1,1] <- base.matrix[r1,1]*x.conc.lowCI[1]/base.matrix[r1,]%*%x.conc.lowCI
	x.lowCI.matrix[r1,2] <- base.matrix[r1,2]*x.conc.lowCI[2]/base.matrix[r1,]%*%x.conc.lowCI

	y.lowCI.matrix[r1,1] <- base.matrix[r1,1]*y.conc.lowCI[1]/base.matrix[r1,]%*%y.conc.lowCI
	y.lowCI.matrix[r1,2] <- base.matrix[r1,2]*y.conc.lowCI[2]/base.matrix[r1,]%*%y.conc.lowCI

	x.hiCI.matrix[r1,1] <- base.matrix[r1,1]*x.conc.hiCI[1]/base.matrix[r1,]%*%x.conc.hiCI
	x.hiCI.matrix[r1,2] <- base.matrix[r1,2]*x.conc.hiCI[2]/base.matrix[r1,]%*%x.conc.hiCI

	y.hiCI.matrix[r1,1] <- base.matrix[r1,1]*y.conc.hiCI[1]/base.matrix[r1,]%*%y.conc.hiCI
	y.hiCI.matrix[r1,2] <- base.matrix[r1,2]*y.conc.hiCI[2]/base.matrix[r1,]%*%y.conc.hiCI
	

	if(r1==1) {
	  xedge.matrix[r1,1] <- base.matrix[r1,1]*x.conc[1]/base.matrix[r1,]%*%x.conc
	  xedge.matrix[r1,2] <- base.matrix[r1,2]*x.conc[2]/base.matrix[r1,]%*%x.conc
	  
	  yedge.matrix[r1,1] <- base.matrix[r1,1]*y.conc[1]/base.matrix[r1,]%*%y.conc
	  yedge.matrix[r1,2] <- base.matrix[r1,2]*y.conc[2]/base.matrix[r1,]%*%y.conc
	  
	}
    }
  }#end r loop

  x.basic <- x.matrix%*%x.points
  y.basic <- y.matrix%*%y.points
  
  x.edge <- xedge.matrix%*%x.points
  y.edge <- yedge.matrix%*%y.points

  x.edge.lowCI <- x.lowCI.matrix%*%x.points
  y.edge.lowCI <- y.lowCI.matrix%*%y.points

  x.edge.hiCI 	<- x.hiCI.matrix%*%x.points
  y.edge.hiCI	<- y.hiCI.matrix%*%x.points
  
  x.outer.loCI <- x.matrix%*%x.lowCI
  x.outer.hiCI <- x.matrix%*%x.hiCI

  y.outer.loCI <- y.matrix%*%y.lowCI
  y.outer.hiCI <- y.matrix%*%y.hiCI
  
    plot(x=x.points, y=y.points, xlim=(range(x.points)+ c(-1,1)*2*x.sd), ylim=range(y.points)+  c(-1,1)*2*y.sd,xlab=expression(paste(delta,C[13])), ylab=expression(paste(delta,N[15])),pch=c(19,19),col=c("grey","grey"))

    points(x.basic,y.basic,type='l',lwd=2,col="grey")  
    
  box(lwd=2)
  if(plot.mix ==TRUE) {
    plotCI(x=mixmu.med[,1], y=mixmu.med[,2], liw=(mixmu.med[,1]-mixmu.loCI[,1]), uiw=(mixmu.hiCI[,1] - mixmu.med[,1]), err="x", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col=c("dimgrey"),pch=19)
    plotCI(x=mixmu.med[,1], y=mixmu.med[,2], liw=(mixmu.med[,2]-mixmu.loCI[,2]), uiw=(mixmu.hiCI[,2] - mixmu.med[,2]), err="y", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col=c("dimgrey"),pch=19)
    points(mixmu.med,lwd=2)
    title('Estimates')
  }
  
	rho.vec <- jags.1$BUGSoutput$mean$rho.source
# 	theta <- seq(0, 2 * pi, length=100)
	for(i in 1:2) {
		
		T.mat <- diag(2)
        T.mat[1,2] <- rho.vec[i]
        T.mat[2,1] <- rho.vec[i]

        sd.vec <- c(cov.source[i], cov.source[i+2])    
		T.mat[1,2] <- rho.vec[i]^2 
		T.mat[2,1] <- rho.vec[i]^2 
		
         if(me.flag) {me.mat <- diag(sigmaz.med[1:2]^2)
            T.mat  <- c(cov.source[i], cov.source[i], cov.source[i+2], cov.source[i+2])*T.mat*sd.vec + me.mat
          } else {
            T.mat  <- c(cov.source[i], cov.source[i], cov.source[i+2], cov.source[i+2])*T.mat*sd.vec
        }

        lines(ellipse(T.mat, centre=c(x.points[i], y.points[i]), level=0.95), lty=2, lwd=2)

	}
	
 ##Plots observed istope values of individuals
    if(plot.ind.flag) {  
        if(me.flag) {
            plotCI(x=X[,1], y=X[,2], uiw=2*sqrt(sigmaz.med[1]), err="x", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col='dimgrey', pch=19)
            plotCI(x=X[,1], y=X[,2], uiw=2*sqrt(sigmaz.med[2]), err="y", xlab=expression(paste(delta,C[13])), add=TRUE, ylab=expression(paste(delta,N[15])), col='dimgrey', pch=19)
        }
        points(X,lwd=2)
        title("Observations")
    }  

}#end biplots



#######################################################
##plots smoothed histograms of individual distributions
#######################################################
curves.plot <- function(jags.1, num.sources, num.chains, color=FALSE,individuals, xlab.vec) {

  output <- jags.1$summary$quantiles
  out.mcmc <- jags.1$mcmc
    
  outnames.temp <- dimnames(jags.1$summary$quantiles)
  out.names <- dimnames(output)[[1]]

  out.perc <- which(outnames.temp[[2]] == "50%")
  pIso.ind.array <- array(NA, c(2,individuals,3,3))
  p.ind.array <- array(NA,c(individuals,3,3))
  p1.popmed <- vector('numeric',num.sources)
  
  mcmc.rbind <- NULL
  for(i in 1:num.chains) {
    mcmc.rbind <- rbind(mcmc.rbind, out.mcmc[[i]])
  }

  par(mfrow=c(num.sources,1))

  for(j in 1:num.sources) {

    p.ind <- which(out.names== paste("p[",1,',',j,']',sep='') )
    p.pop <- which(out.names== paste("p.pop[",j,']',sep='') )

    p1.popmed[j] <- output[p.pop, out.perc]

    temp.pop <- hist(mcmc.rbind[,p.pop],plot=FALSE, breaks=seq(0,1,length.out=100))
    pop.smooth <- lowess(temp.pop$mids, y= temp.pop$density,f=0.1)
    
    plot(pop.smooth, type='l', lwd=2, col="black", xlim=c(-0.1,1),ylim=c(0,max(c(pop.smooth$y))*1.1),ylab="probability density", xlab=xlab.vec[j])

    for(i in 2:individuals) {
        
     p.ind <- which(out.names== paste("p[",i,',',j,']',sep='') )
     temp <- hist(mcmc.rbind[,p.ind],plot=FALSE, breaks=seq(0,1,length.out=100))
     temp.smooth <- lowess(c(0,temp$mids,1), y= c(0,temp$density,0), f=.1)
     temp.smooth$x <- temp.smooth$x[-c(1,length(temp.smooth$x))]
     temp.smooth$y <- temp.smooth$y[-c(1,length(temp.smooth$y))]
     
     
     if(color) { lines(temp.smooth, type='l', lwd=2, col=rgb(0,0,128,alpha=100,maxColorValue=255)) } else {
     lines(temp.smooth, type='l', lwd=1, col="white")}    
  }
  
  if(min(pop.smooth$x) > 0.025) {pop.smooth$x <- c(0.025,pop.smooth$x); pop.smooth$y <- c(0,pop.smooth$y)}
  if(max(pop.smooth$x) < 0.975) {pop.smooth$x <- c(pop.smooth$x,0.975); pop.smooth$y <- c(pop.smooth$y,0)}
  lines(pop.smooth, type='l', lwd=2, col="black")
   
  }#end for
}#end curves.plot
