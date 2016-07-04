\name{flomax}
\alias{flomax}
\title{

  ML estimation of the Lomax distribution

}
\description{
  
  Fast Maximum Likelihood estimation of the Lomax distribution.

}
\usage{
flomax(x,
       info.observed = TRUE,
       plot = FALSE,
       scaleData = TRUE,
       cov = TRUE)
}
\arguments{
  \item{x}{

    Sample vector to be fitted. Should contain only positive
    non-NA values.

  }
  \item{info.observed}{

    Should the observed information matrix be used or the expected one
    be used?

  }
  \item{plot}{

    Logical. If \code{TRUE}, a plot will be produced showing the
    derivative of the concentrated log-likelihood, function of the shape
    parameter.

  }
  \item{scaleData}{

    Logical. If \code{TRUE} observations in \code{x} (which are
    positive) are divided by their mean value.  The results are in
    theory not affected by this transformation, but scaling the data
    could improve the estimation in some cases.  The log-likelihood
    plots are shown using the scaled values so the returned estimate of
    the scale parameter is not the the abscissa of the maximum shown on
    the plot.

  }
  \item{cov}{
    
     Logical. If \code{FALSE}, a minimal estimation is performed with no
     covariance matrix or derivative returned. This can be useful when a
     large number of ML estimations are required, e.g. to sample from a
     likelihood ratio.

   }
}
\details{

  The likelihood is concentrated with respect to the shape parameter.
  This function is increasing for small values of the scale parameter
  \eqn{\beta}. For large \eqn{\beta}, the derivative of the concentrated
  log-likelihood tends to zero, and its sign is that of \eqn{(1 -
  \textrm{CV}^2)}{(1-CV^2)} where \eqn{\textrm{CV}}{CV} is the
  coefficient of variation, computed using \eqn{n} as denominator in the
  formula for the standard deviation.
  
  The ML estimate does not exist when the sample has a coefficient of
  variation CV less than 1 and it may fail to be found when CV is
  greater than yet close to 1.
  
}
\value{
   A list with the following elements
   \item{estimate}{

     Parameter ML estimates.

   }
   \item{sd}{

     Vector of (asymptotic) standard deviations for the estimates.

   }
   \item{loglik}{

     The maximised log likelihood.

   }
   \item{dloglik}{

     Gradient of the log-likelihood at the optimum. Its two
     elements should normally be close to zero.

   }
   \item{cov}{

     The (asymptotic) covariance matrix computed from theoretical or
     observed information matrix.

   }
   \item{info}{

     The information matrix.

   }
}
\references{
  J. del Castillo and J. Daoudi (2009) "Estimation of the Generalized
  Pareto Distribution", \emph{Statist. Probab. Lett.} 79(5),
  pp. 684-688.

  D.E. Giles, H. Feng & R.T. Godwin (2013) "On the Bias of the Maximum
  Likelihood Estimator for the Two-Parameter Lomax Distribution"
  \emph{Comm. Statist. Theory Methods}. Vol. 42, n. 11, pp. 1934-1950.

  N. Johnson, S. Kotz and N. Balakrishnan 
  \emph{Continuous Univariate Distributions} vol. 1, Wiley 1994.

}
\author{

  Yves Deville
}
\note{

  The estimates are biased for small or medium sized sample. The bias is
  positive for the shape parameter, thus the estimated shape tends to be
  larger than the true unknown value.

  Fitting a Lomax distribution to an exponential sample might lead to a
  divergence since the exponential is the limit of a Lomax distribution
  with large shape and large scale with constant ratio shape/scale.
  Fitting this distribution to a sample having a coefficient of
  variation smaller than 1 is not allowed since it should lead to
  divergence of the estimation.

  The default value of \code{info.observed} was set to \code{TRUE} from
  version \code{3.0-1} because standard deviations obtained with this
  choice are usually better.
  
}
\seealso{
  
  \code{\link{Lomax}} for the Lomax distribution.

}
\examples{
## generate sample
set.seed(1234)
n <- 200
alpha <- 2 + rexp(1)
beta <- 1 + rexp(1)
x <- rlomax(n, scale = beta, shape = alpha)
res <- flomax(x, plot = TRUE)

## compare with a GPD with shape 'xi' and scale 'sigma'
xi <- 1 / alpha; sigma <- beta * xi  
res.evd <- evd::fpot(x, threshold = 0, model = "gpd")
xi.evd <- res.evd$estimate["shape"]
sigma.evd <- res.evd$estimate["scale"]
beta.evd <- sigma.evd / xi.evd 
alpha.evd <- 1 / xi.evd
cbind(Renext = res$estimate, evd = c(alpha = alpha.evd, beta = beta.evd))
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
