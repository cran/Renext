\name{ini.mixexp2}
\alias{ini.mixexp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Simple  estimation for the mixture of two exponential distributions

}
\description{

  Compute a simple (preliminary) estimation for the tree parameters of
  the mixture of two exponential distributions

}
\usage{
   ini.mixexp2(x, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{

    Sample: numerical vector with elements \code{>0}.

  }
  %%  \item{pct}{Pair of percents for the estimation. The upper percent
  %%    indicates the largest data to estimate \code{prob1} and \code{rate1}
  %%    by linear regression. The lower indicates the smallest data to
  %%    estimate \code{rate2} using a regression line trhough the origin.}
  \item{plot}{
    
    Should a graphic be displayed?

  }
}
\details{
  
  This function gives estimators using several methods if necessary. The
  goal is to find the rates \code{rate1}, \code{rate2} and the mixing
  probability \code{prob1} with the 'feasibility' constraints \code{0 <
  rate1} \code{< rate2} and \code{0 < prob1 < 1}.

  First the method of moments is used. If the estimates are feasible
  they are returned with \code{method} \code{=} \code{"moments"}.  If
  not, the estimates are derived using two linear regressions. A
  regression without constant using only the smallest values gives an
  estimator of the mean rate.  A regression using only the largest
  values gives \code{rate1} and \code{prob1}. Yet the constraints must
  be fulfilled. If they are, the estimates are returned (together with
  \code{method =} \code{"Hreg"} suggesting a cumulative hazard
  regression). If not, a (poor) default estimate is returned with
  \code{method =} \code{"arbitrary"}.
  
}
\value{
  A list 
  \item{estimate}{

    A vector with named elements \code{"prob1"}, \code{"rate1"} and
    \code{"rate2"}.

  }
  \item{method}{

    The method that really produced the estimators.

  }
}
%%\references{ }
\author{Yves Deville}
\note{  

  The method of moments is implemented in \code{mom.mixexp2}.  Further
  investigations are needed to compare the estimators (moments or Hreg)
  and select the best strategy.

  Note that this function returns the estimate within a list and no
  longer as a vector with named elements as was the case before.

}
\seealso{See \code{\link{MixExp2}}, \code{\link{mom.mixexp2}}. }
\examples{
set.seed(1234)
x <- rmixexp2(n = 100, prob1 = 0.5, rate2 = 4)
res <- ini.mixexp2(x, plot = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
