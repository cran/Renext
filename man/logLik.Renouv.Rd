\name{logLik.Renouv}
\alias{AIC.Renouv}
\alias{BIC.Renouv}
\alias{logLik.Renouv}
\alias{nobs.Renouv}

\title{
  Log-likelihood of a "Renouv" object
}
\description{
  Log-likelihood, AIC, BIC and number of observations of an object of class "Renouv".
}
\usage{
\S3method{AIC}{Renouv}(object, ..., k = 2)
\S3method{BIC}{Renouv}(object, ...)
\S3method{logLik}{Renouv}(object, ...)
\S3method{nobs}{Renouv}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    
    Object of class \code{"Renouv"}.

  }
  \item{k}{

    See \code{\link{AIC}}.

  }
  \item{\dots}{

    Not used yet.

  }
}
%% \details{}
%% \value{}
%%\references{ }
\author{
  Yves Deville
}
\note{

  \code{logLik}, \code{AIC} and \code{BIC} can be used with an object of
  class \code{"Renouv"} which makes use of historical data. In this
  case, the number of observations may be misleading since a single
  historical observation may concern dozens of years and thus have a
  much greater impact on the estimation of the tail than an "ordinary"
  observation.

}
\section{Caution}{
  
  Comparing log-likelihoods, AIC or BIC for different \code{Renouv}
  objects makes sense only when these share the same data and the same
  threshold.

} 
\seealso{
  The \code{\link{AIC}}, \code{\link{nobs}} generic functions.
}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
