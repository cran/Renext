\name{GPD}
\alias{GPD}
\alias{dGPD}
\alias{pGPD}
\alias{qGPD}
\alias{rGPD}
\alias{hGPD}
\alias{HGPD}
\title{Generalised Pareto Distribution}
\description{
  
  Density function, distribution function, quantile function, random
  generation, hazard and cumulative hazard functions for the Generalised
  Pareto Distribution.

}
\usage{
   dGPD(x, loc = 0.0, scale = 1.0, shape = 0.0, log = FALSE)
   pGPD(q, loc = 0.0, scale = 1.0, shape = 0.0, lower.tail = TRUE)
   qGPD(p, loc = 0.0, scale = 1.0, shape = 0.0, lower.tail = TRUE)
   rGPD(n, loc = 0.0, scale = 1.0, shape = 0.0)
   hGPD(x, loc = 0.0, scale = 1.0, shape = 0.0)
   HGPD(x, loc = 0.0, scale = 1.0, shape = 0.0) 
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
  \item{loc}{

    Location parameter \eqn{\mu}.

  }
 \item{scale}{

    Scale parameter \eqn{\sigma}.

  }
  \item{shape}{

    Shape parameter \eqn{\xi}.

  }
  \item{log}{

    Logical; if \code{TRUE}, the log density is returned.

  }
  \item{lower.tail}{

    Logical; if \code{TRUE} (default), probabilities are
    \eqn{\textrm{Pr}[X <= x]}{Pr[X <= x]}, otherwise, \eqn{\textrm{Pr}[X
    > x]}{Pr[X > x]}.

  }
}
\details{
  Let \eqn{\mu}, \eqn{\sigma} and \eqn{\xi} denote \code{loc},
  \code{scale} and \code{shape}. The distribution values \eqn{y}
  are \eqn{\mu \leq y < y_{\textrm{max}}}{\mu <= y < ymax}.
  
  When \eqn{\xi \neq 0}{\xi != 0}, the survival function value for
  \eqn{y \geq \mu}{y >= \mu} is given by

  \deqn{S(y) = \left[1 + \xi(y - \mu)/\sigma\right]^{-1/ \xi} \qquad
    \mu < y < y_{\textrm{max}}}{ S(y) = [1 + \xi (y - \mu) /
    \sigma]^(-1 / \xi) (\mu < y < ymax)}
  where the upper end-point is
    \eqn{y_{\textrm{max}} = \infty}{ymax = Inf} for
    \eqn{\xi >0} and \eqn{y_{\textrm{max}} = \mu -\sigma/ \xi}{ymax =
    \mu -\sigma / \xi} for \eqn{\xi <0}.

  When \eqn{\xi = 0}, the distribution is exponential with survival
   \deqn{S(y) = \exp\left[- (y - \mu)/\sigma\right] \qquad \mu \leq y. }{
    S(y) = exp(-(x-\mu) / \sigma)      (\mu <= y)}
  
}
\value{
  
  \code{dGPD} gives the density function, \code{pGPD} gives the
  distribution function, \code{qGPD} gives the quantile function, and
  \code{rGPD} generates random deviates. The functions
  \code{hGPD} and \code{HGPD} return the hazard rate and the cumulative
  hazard.
  
}
\note{

  The functions are slight adaptations of the \code{[r,d,p,q]gpd}
  functions in the \pkg{evd} package. The main difference is that
  these functions return \code{NaN} when \code{shape} is negative, as
  it might be needed in unconstrained optimisation. The quantile function
  can be used with \code{p=0} and \code{p=1}, then returning the lower and
  upper end-point.       
}
%%\references{}

%%\references{  }	
\seealso{

  \code{\link{fGPD}} to fit such a distribution by Maximum Likelihood.

}
\examples{
qGPD(p = c(0, 1), shape = -0.2)
shape <- -0.3
xlim <- qGPD(p = c(0, 1), shape = shape)
x <- seq(from = xlim[1], to = xlim[2], length.out = 100)
h <- hGPD(x, shape = shape)
plot(x, h, type = "o", main = "hazard rate for shape < 0")
shape <- 0.2
xlim <- qGPD(p = c(0, 1 - 1e-5), shape = shape)
x <- seq(from = xlim[1], to = xlim[2], length.out = 100)
h <- hGPD(x, shape = shape)
plot(x, h, type = "o", main = "hazard rate shape > 0 ")
}
\keyword{distribution}

