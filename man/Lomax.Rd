\name{Lomax}
\alias{Lomax}
\alias{dlomax}
\alias{plomax}
\alias{qlomax}
\alias{rlomax}
\title{Lomax distribution}
\description{
  
  Density function, distribution function, quantile function and
  random generation for the Lomax distribution. 

}
\usage{
   dlomax(x, scale = 1.0, shape = 4.0, log = FALSE)
   plomax(q, scale = 1.0, shape = 4.0, lower.tail = TRUE)
   qlomax(p, scale = 1.0, shape = 4.0)
   rlomax(n, scale = 1.0, shape = 4.0) 
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
  \item{scale, shape}{

    Scale and shape parameters. Vectors of length > 1 are not accepted.

  }
  \item{log}{

    Logical; if \code{TRUE}, the log density is returned.

  }
  \item{lower.tail}{

    Logical; if \code{TRUE} (default), probabilities are
    \eqn{\textrm{Pr}[X <= x]}{Pr[X <= x]}, otherwise, \eqn{\textrm{Pr}[X
    > x]}{Pr[X > x]}.}

}
\details{
  
  The Lomax distribution function with shape \eqn{\alpha > 0} and scale
  \eqn{\beta > 0} has survival function
  \deqn{S(y) = \left[1 + y/\beta \right]^{-\alpha} \qquad (y > 0)}{
    S(y) = [1 + y / \beta ]^(-\alpha) (y > 0)}
  This distribution has increasing hazard and decreasing mean
  residual life (MRL). The coefficient of variation decreases with
  \eqn{\alpha}, and tends to \eqn{1}{1} for large \eqn{\alpha}.  The
  default value \eqn{\alpha=4} corresponds to \eqn{\textrm{CV} =
  \sqrt{2}}{CV = sqrt(2)}.

}
\value{

  \code{dlomax} gives the density function, \code{plomax} gives the
  distribution function, \code{qlomax} gives the quantile function, and
  \code{rlomax} generates random deviates.

}
\note{

  This distribution is sometimes called \emph{log-exponential}. It is a
  special case of Generalised Pareto Distribution (GPD) with positive
  shape \eqn{\xi > 0}, scale \eqn{\sigma} and location \eqn{\mu=0}.  The
  Lomax and GPD parameters are related according to \deqn{\alpha =
  1/\xi, \qquad \beta = \sigma/\xi.}{\alpha = 1/\xi, \beta =
  \sigma/\xi.}  The Lomax distribution can be used in POT to describe
  excesses following GPD with shape \eqn{\xi>0} thus with decreasing
  hazard and increasing Mean Residual Life.
  
  Note that the exponential distribution with rate \eqn{\nu} is the
  limit of a Lomax distribution having large scale \eqn{\beta} and large
  shape \eqn{\alpha}, with the constraint on the shape/scale ratio
  \eqn{\alpha/\beta = \nu}.
  
}
\references{
  
  Johnson N. Kotz S. and N. Balakrishnan 
  \emph{Continuous Univariate Distributions} vol. 1, Wiley 1994.

  \href{https://en.wikipedia.org/wiki/Lomax_distribution}{Lomax distribution in Wikipedia}
  
}

%%\references{  }	
\seealso{

  \code{\link{flomax}} to fit the Lomax distribution by Maximum
  Likelihood. 

}

\examples{
shape <- 5; scale <- 10
xl <- qlomax(c(0.00, 0.99), scale = scale, shape = shape)
x <- seq(from = xl[1], to = xl[2], length.out = 200)
f <- dlomax(x, scale = scale, shape = shape)
plot(x, f, type = "l", main = "Lomax density")
F <- plomax(x, scale = scale, shape = shape)
plot(x, F, type ="l", main ="Lomax distribution function")
}
\keyword{distribution}

