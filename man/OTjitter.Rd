\name{OTjitter}
\alias{OTjitter}
\title{
  Add a small amount of noise to a numeric vector
}
\description{
  Add a small amount of noise to a numeric vector keeping all the
  values above the given threshold.
}
\usage{
OTjitter(x, threshold = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The numeric vector to which \emph{jitter} should be added.
  }
  \item{threshold}{A threshold above which all elements of the
    modified vector must stay.
  }
}
%%\details{}
\value{
  A vector with the same length and nearly the same values
  as \code{x}. As in \code{\link{jitter}}, a small amount of noise is
  added to each value of \code{x}. The noise level is adjusted so that
  every noisy value remains above the specified threshold. When the
  a value is very close to the threshold, only a very small amount of
  negative noise can be added.
}
\note{
  The aim of this function is to remove possible ties in experimental OT
  data. Ties cause problems or warnings in some goodness-of-fit tests such
  as Kolmogorov-Smirnov.
}
%%\references{}
\author{Yves Deville}
\seealso{
  \code{\link{jitter}}
}
\examples{
## Garonne data (heavily rounded)
x <- Garonne$OTdata$Flow
min(x) 
xmod <- OTjitter(x, threshold = 2500)
length(x)
nlevels(as.factor(x))
nlevels(as.factor(xmod))
max(abs(x-xmod))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
