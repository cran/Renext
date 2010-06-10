\name{Brest.years.missing}
\alias{Brest.years.missing}
\docType{data}
\title{Years with missing periods in 'Brest.year' dataset}
\description{
  Years with missing periods in the 'Brest.years' dataset
}
\usage{data(Brest.years.missing)}
\format{
  The format is:
 int [1:49] 1846 1847 1852 1857 1858 1859 1860 1861 1862 1863 ...
}
\details{
  Vector of years containing missing periods in the
  \code{\link{Brest.years}} dataset. This years should be ignored
  when computing yearly statistics sucha as event rates since sime
  records are lost.
}
%%\source{}
%%\references{}
\examples{
data(Brest.years.missing)
print(Brest.years.missing)
}
\keyword{datasets}
