\name{SLTW}
\alias{SLTW}
\alias{dSLTW}
\alias{pSLTW}
\alias{qSLTW}
\alias{rSLTW}
\title{Shifted Left Truncated Weibull (SLTW) distribution}
\description{
  
  Density function, distribution function, quantile function and random
  generation for the Shifted Left Truncated Weibull distribution.

}
\usage{
   dSLTW(x, delta = 1.0, shape = 1.0, scale = 1.0, log = FALSE)
   pSLTW(q, delta = 1.0, shape = 1.0, scale = 1.0, lower.tail = FALSE)
   qSLTW(p, delta = 1.0, shape = 1.0, scale = 1.0)
   rSLTW(n, delta = 1.0, shape = 1.0, scale = 1.0) 
}
\arguments{
  \item{x, q}{

    Vector of quantiles.

  }
  \item{p}{

    Vector of probabilities.

  }
  \item{n}{
    Number of observations.

  }
  \item{delta, shape, scale}{

    Shift, shape and scale parameters. Vectors of length > 1 are not
    accepted.

  }
  \item{log}{

    Logical; if \code{TRUE}, the log density is returned.

  }
  \item{lower.tail}{

    Logical; if \code{TRUE} (default), probabilities are
    \eqn{\textrm{Pr}[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{\textrm{Pr}[X > x]}{Pr[X > x]}.

  }
}
\details{
  
  The SLTW distribution function with shape \eqn{\alpha > 0}, scale
  \eqn{\beta > 0} and shift \eqn{\delta > 0} has survival function
  \deqn{S(y) = 
    \exp\left\{
       -\left[  \left( \frac{y + \delta}{\beta} \right)^\alpha
       -        \left( \frac{\delta}{\beta}   \right)^\alpha
       \right]
    \right\} \qquad (y > 0)
  }{
    S(y) = exp{- [(y + \delta)/\beta]^\alpha + [\delta/\beta]^\alpha}
              (y > 0)
  }
  This distribution is that of \eqn{Y := X - \delta} conditional to
  \eqn{X > \delta} where \eqn{X} follows a Weibull distribution with
  shape \eqn{\alpha} and scale \eqn{\beta}.
  
  The hazard and mean residual life (MRL) are monotonous functions with
  the same monotonicity as their Weibull equivalent (with the same shape
  and scale). The moments or even expectation do not have simple
  expression.

  This distribution is sometimes called \emph{power exponential}. It is
  occasionally used in POT with the shift \code{delta} taken as the
  threshold as it should be when the distribution for the level \eqn{X}
  (and not for the exceedance \eqn{Y}) is known to be the standard
  Weibull distribution.  } \value{ \code{dSLTW} gives the density
  function, \code{pSLTW} gives the distribution function, \code{qSLTW}
  gives the quantile function, and \code{rSLTW} generates random
  deviates.

}
%%\references{  }	
\seealso{

  \code{\link{Lomax}} for the Lomax distribution which is a limit case
  of \code{SLTW}.

}

\examples{
shape <- rexp(1)+1  
delta = 10

xl <- qSLTW(c(0.001, 0.99), delta = delta, shape = shape)
x <- seq(from = xl[1], to = xl[2], length.out = 200)
f <- dSLTW(x, delta = delta, shape = shape)
plot(x, f, type = "l", main = "SLTW density")

F <- pSLTW(x, delta = delta, shape = shape)
plot(x, F, type = "l", main = "SLTW distribution")

X <- rSLTW(5000, delta = delta, shape = shape)
## Should be close to the uniform repartition
plot(ecdf(pSLTW(X, delta = delta, shape = shape)))

}
\keyword{distribution}

