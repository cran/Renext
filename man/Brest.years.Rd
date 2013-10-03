\name{Brest.years}
\alias{Brest.years}
\docType{data}
\title{Surge heights at Brest partial data}
\description{
  Surge heights at Brest (France) 
}
\usage{Brest.years}
\format{
  A data frame with 954 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{Year e.g; 1980}
    \item{\code{Surge}}{Surge heights above the threshold of 30 cm.}
  }
}
\details{
  These data are a simplified version of \code{\link{Brest}}. For each
  surge event only the year is retained as timestamp. Years with
  missing periods are available as a vector \code{\link{Brest.years.missing}}. 
  
  This dataset is useful for testing since similar data are sometimes
  met in the analyses.
  
}
%\source{ }
%%\references{}
\examples{
names(Brest.years)
}
\keyword{datasets}
