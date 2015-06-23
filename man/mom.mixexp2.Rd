\name{mom.mixexp2}
\alias{mom.mixexp2}
\title{

  Moment estimation for the mixture of two exponential distributions

}
\description{
  
  Compute the moment estimation for the tree parameters of the mixture of
  two exponential distributions

}
\usage{
   mom.mixexp2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{

    Sample. Vector containing values \code{>0}.

  }
}
\details{
  
  The three parameters (probability and the two rates) are computed from
  the first three moments (theoretical and sample). It can be shown that
  the inverse rates are obtained solving a quadratic equation. However
  the roots can be negative or complex and the estimates are not valid
  ones.

}
\value{
  A list with elements
  \item{estimate}{
    
    A vector with named elements \code{"prob1"}, \code{"rate1"} and
    \code{"rate2"}. When the moment estimators are not valid (negative
    or complex rates), a vector of three \code{NA} is returned.

  }
  \item{method}{
    
    Character \code{"moments"}.

  }  
}
\references{

  Paul R. Rider. The Method of Moments Applied to a Mixture of Two
  Exponential Distributions. \emph{Ann. Math. Statist.} Vol. 32, Number
  1 (1961), 143-147.

}
\author{

  Yves Deville

}
\note{

  The theoretical coefficient of variation (CV) of a mixture of two
  exponential distributions always exceeds 100\%. When the sample CV is
  smallest than 100\%, no valid estimates exist since the two first
  moments can not be matched.

}
\seealso{

  See \code{\link{ini.mixexp2}} for a more versatile initial estimation.

}
\examples{

x <- rmixexp2(n = 100, prob1 = 0.5, rate1 = 1.0, rate2 = 3.0)
est <- mom.mixexp2(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
