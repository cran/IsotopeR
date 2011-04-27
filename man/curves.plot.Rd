\name{curves.plot}
\alias{curves.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Curves Plot
}
\description{
Plots the population level proportional source estimates (black) and individual level estimates (blue).
}
\usage{
curves.plot(jags.1, num.sources, num.chains, color = FALSE, individuals, xlab.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jags.1}{
A jags object, output from the runjags package.
}
  \item{num.sources}{
Number of isotopic sources used
}
  \item{num.chains}{
Number of mcmc chains
}
  \item{color}{
Whether to make plots in color
}
  \item{individuals}{
Number of individuals in the sample
}
  \item{xlab.vec}{
The name of dietary sources.
}
}

\author{
Jake Ferguson
}


\examples{
#IsotopeR()
#load(SampleOutput.Rdata)
#curves.plot(jags.1=jags.out, num.sources=3, individuals=7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
