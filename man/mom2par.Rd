\name{mom2par}
\alias{mom2par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameters from moments}
\description{
  Compute parameters from (theoretical) moments
}
\usage{
   mom2par(densfun = "exponential",
           mean,
           sd = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{densfun}{Name of the distribution. This can be at present time:
    \code{"exponential"}, \code{"weibull"}, \code{"gpd"},
    \code{"gamma"}, \code{"negative binomial"}.
  }
  \item{mean}{Theoretical mean (expectation) of the distribution. Can be
    a vector, in which case the parameters will be vectors.}
  \item{sd}{Standard deviation.}
}
\details{
  For some distributions like Weibull, it is necessary to find
  a numerical solution since the parameters have no closed form
  expression involving the moments.
}
\value{
  A named list containing the parameters values e.g. with names \code{shape}
  and \code{scale}. When \code{mean} or \code{sd} is vector the list
  contains vectors. 
}
%%\references{}
\author{Yves Deville}
\note{
  The name of the formal argument \code{densfun} is for compatibility with
  \code{fitdistr} from the MASS package. However, unlike in
  \code{fitdistr}  this formal can not be given a density value, i.e. an
  object of the class \code{"function"} such as \code{dnorm}. 
}
%%\seealso{}
\examples{
## Weibull
mom2par(densfun = "weibull", mean = 1, sd = 2) 
## Genrealised Pareto
mom2par(densfun = "gpd", mean = 1, sd = 2)
## Gamma
mom2par(densfun = "gamma", mean = 1:10, sd = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
