\name{Tri.plots}
\alias{Tri.plots}
\title{
Mixing Space, 3 sources
}
\description{
3 source mixing space plot used to plot the source estimations and the observations or model fits.
}
\usage{
Tri.plots(jags.1, X, sources=NA, plot.mix = FALSE, plot.ind.flag = FALSE, me.flag = FALSE, color.plots=TRUE)
}
\arguments{
\item{jags.1}{A jags object, output from the runjags package.}
\item{X}{Isotopic mixture data}
\item{sources}{A matrix including all source information. Should be formatted as a source file as discussed in the IsotopeR vignette.}
\item{plot.mix}{A flag, plots the estimated mixture values, when TRUE (1). If FALSE (0) then the the function plots the observed mixture values. }
\item{plot.ind.flag}{Legacy argument, not currently used.}
\item{me.flag}{Denotes wether to plot the measurement error along with observations.}
\item{color.plots}{ Denotes whether plots should be color, or greyscale. The user does not have the ability to specify which colors should be used. }
}
\author{
Jake Ferguson
}
\examples{
#IsotopeR()
#load(SampleOutput.Rdata)
#Tri.plots(jags.1=jags.out, X=X)
}

