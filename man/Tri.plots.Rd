\name{Tri.plots}
\alias{Tri.plots}
\title{
Mixing Space, 3 sources
}
\description{
3 source mixing space plot used to plot the source estimations and the observations or model fits.
}
\usage{
Tri.plots(jags.1, X, plot.mix = FALSE, plot.ind.flag = FALSE, me.flag = FALSE)
}
\arguments{
  \item{jags.1}{A jags object, output from the runjags package.}
  \item{X}{Isotopic mixture data}
  \item{plot.mix}{A flag, plots the estimated mixture values, when TRUE (1). If FALSE (0) then the the function plots the observed mixture values. }
  \item{plot.ind.flag}{Legacy argument, not currently used.}
  \item{me.flag}{Denotes wether to plot the measurement error along with observations.}
}
\author{
Jake Ferguson
}
\examples{
#IsotopeR()
#load(SampleOutput.Rdata)
#Tri.plots(jags.1=jags.out, X=X)
}

