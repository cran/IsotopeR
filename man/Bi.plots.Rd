\name{Bi.plots}
\alias{Bi.plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixing Space, 2 sources
}
\description{
2 source mixing space plot used to plot the source estimations and the observations or model fits.
}
\usage{
Bi.plots(jags.1, X, plot.mix = FALSE, plot.ind.flag = FALSE, me.flag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jags.1}{
A jags object, output from the runjags package.
}
  \item{X}{
Isotopic mixture data
}
  \item{plot.mix}{
A flag, plots the estimated mixture values, when TRUE (1). If FALSE (0) then the the function plots the observed mixture values. 
}
  \item{plot.ind.flag}{
Legacy argument, not currently used.
}
  \item{me.flag}{
Denotes wether to plot the measurement error along with observations.
}
}
\author{
Jake Ferguson
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
#IsotopeR()
#load(SampleOutput.Rdata)
#Bi.plots(jags.1=jags.out, X=X)
}

