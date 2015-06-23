\name{vcov.Renouv}
\alias{vcov.Renouv}

\title{
  
  Variance-covariance matrix of the estimates of a "Renouv" object

}
\description{

  Variance-covariance matrix of the estimates of a "Renouv" object.

}
\usage{
\S3method{vcov}{Renouv}(object, ...)
}
\arguments{
  \item{object}{

    Object of class \code{"Renouv"}.

  }
  \item{\dots}{

    Not used at the time.

  }
}
%% \details{}
\value{
  
  A variance-covariance matrix. The rows an columns correspond to the
  parameters of the Renouv object. The are the rate \code{"lambda"} for
  the Poisson process, and the parameters of the distribution for the
  excesses over the threshold, with names depending on the chosen
  distribution.
  
}
\author{

  Yves Deville

}
\seealso{

  The \code{\link{vcov}} generic.

}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
