\name{fmaxlo}
\alias{fmaxlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ML estimation of a 'maxlo' distribution}
\description{
  Fast Maximum Likelihood estimation of a 'maxlo' distribution
}
\usage{
fmaxlo(x,
       info.observed = FALSE,
       plot = FALSE,
       scaleData = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Sample vector to be fitted. Should contain only positive
    non-NA values.
  }
  \item{info.observed}{Should the observed information matrix be used or
    the expected one be used?
  }
  \item{plot}{Logical. If \code{TRUE}, a plot will be produced
    showing the derivative of the concentrated log-likelihood,
    function of the shape parameter.
  }
  \item{scaleData}{Logical. If \code{TRUE} observations in \code{x}
    (which are positive) are divided by their mean value.
    The results are in theory not affected by this transformation, but
    scaling the data could improve the estimation in some cases. 
    The log-likelihood plots are shown using the scaled values so
    the returned estimate of the scale parameter is not the
    the abscissa of the maximum shown on the plot.
  }
}
\details{
  The 'maxlo' likelihood is concentrated with respect to the shape
  parameter, thus the function to be maximised has only one one scalar argument:
  the scale parameter \eqn{\beta}. For large scale \eqn{\beta},
  the derivative of the concentrated log-likelihood tends
  to zero, and its sign is that of \eqn{(\textrm{CV}^2-1)}{(CV^2-1)} where
  \eqn{\textrm{CV}}{CV} is the coefficient of variation, computed using
  \eqn{n} as denominator in the formula for the standard deviation.
  
  The ML estimate does  not exist when the sample has a coefficient of
  variation CV greater than 1 and it may fail to be found when CV is
  smaller than yet close to 1.
  
  The information matrix can be computed by noticing that when the
  r.v. \eqn{Y}{Y} follows the 'maxlo' distribution with shape
  \eqn{\alpha} and scale \eqn{\beta} the r.v
  \eqn{V:= 1/(1-Y/\beta)} follows a Pareto
  distribution with minimum 1 and and shape parameter
  \eqn{\alpha}{\alpha}. The information matrix involves the second
  order moment of \eqn{V}{V}. 
}
\value{
  A list with the following elements
  \item{estimate}{Parameter ML estimates.}
  \item{sd}{Vector of (asymptotic) standard deviations for the estimates.}
  \item{loglik}{The maximised log-likelihood.}
  \item{dloglik}{Gradient of the log-likelihood at the optimum. Its two
    elements should normally be close to zero.}
  \item{cov}{The (asymptotic) covariance matrix computed from theoretical or
    observed information matrix.}
  \item{info}{The information matrix.}
}
%%\references{}
\author{
  Yves Deville
}
\note{
  The name of the distribution hence also that of the fitting function
  are still experimental and might be changed.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  \code{\link{Maxlo}} for the description of the distribution.
}
\examples{
## generate sample
set.seed(1234)
n <- 200
alpha <- 2 + rexp(1)
beta <- 1 + rexp(1)
x <- rmaxlo(n, scale = beta, shape = alpha)
res <- fmaxlo(x, plot = TRUE)

## compare with a GPD with shape 'xi' and scale 'sigma'
xi <- -1 / alpha; sigma <- -beta * xi
res.evd <- evd::fpot(x, threshold = 0, model = "gpd")
xi.evd <- res.evd$estimate["shape"]
sigma.evd <- res.evd$estimate["scale"]
beta.evd <- -sigma.evd / xi.evd 
alpha.evd <- -1 / xi.evd
cbind(Renext = res$estimate, evd = c(alpha = alpha.evd, beta = beta.evd))
  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
