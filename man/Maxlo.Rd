\name{Maxlo}
\alias{Maxlo}
\alias{dmaxlo}
\alias{pmaxlo}
\alias{qmaxlo}
\alias{rmaxlo}
\title{'maxlo' distribution}
\description{
  Density function, distribution function, quantile function and
  random generation for the 'maxlo' distribution. 
}
\usage{
   dmaxlo(x, scale = 1.0, shape = 4.0, log = FALSE)
   pmaxlo(q, scale = 1.0, shape = 4.0, lower.tail = TRUE)
   qmaxlo(p, scale = 1.0, shape = 4.0)
   rmaxlo(n, scale = 1.0, shape = 4.0) 
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

    Shift and shape parameters. Vectors of length
    > 1 are not accepted.

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
  
  The 'maxlo' distribution function with shape \eqn{\alpha>0} and scale
  \eqn{\beta>0} is a special case of Generalised Pareto (GPD) with
  \emph{negative shape} \eqn{\xi < 0} and location at zero. This is the
  finite upper endpoint case of the GPD. Its name is nonstandard and was
  chosen to suggest some form of symmetry with respect to the Lomax
  distribution.
  
  The survival  function is
  \deqn{S(y) = \left[1-y/\beta\right]^\alpha \qquad 0 < y < \beta}{
    S(y) = [1 - y / \beta]^\alpha      (0 < y < \beta)}

  This distribution has a coefficient of variation smaller than \eqn{1}{1}.
}
\value{
  
  \code{dmaxlo} gives the density function, \code{pmaxlo} gives the
  distribution function, \code{qmaxlo} gives the quantile function, and
  \code{rmaxlo} generates random deviates.

}
\note{

  The 'maxlo' and GPD parameters are related according to \deqn{\alpha =
  -1/\xi, \qquad \beta = -\sigma/\xi.}{\alpha = -1/\xi, \beta =
  -\sigma/\xi.}  where \eqn{\sigma} is the scale parameter of the
  GPD. Since only GPD with \eqn{\xi > -0.5} seem to be used in practice,
  this distribution should be used with \eqn{\alpha > 2}.
  
  This distribution can be used in POT to describe bounded excesses
  following GPD with shape \eqn{\xi < 0}. The scale parameter
  \eqn{\beta} then represents the upper end-point of the excesses,
  implying the finite upper end-point \eqn{u + \beta} for the levels,
  where \eqn{u} is the threshold. It can be used in \code{\link{Renouv}}
  with a fixed scale parameter, thus allowing a control of the upper
  end-point.
   
  This distribution is simply a rescaled version of a beta distribution
  and also a rescaled version of a Kumaraswamy distribution. The name
  "maxlo" is used here to suggest a form of symmetry to Lomax
  distribution.
       
}
%%\references{}

%%\references{  }	
\seealso{

  \code{\link{fmaxlo}} to fit such a distribution by Maximum Likelihood.

}
\examples{
xs <- rmaxlo(500, shape = 2.2, scale = 1000)
hist(xs, main = "'maxlo' distribution"); rug(xs)

xs <- rmaxlo(500, shape = 4, scale = 1000)
hist(xs, main = "'maxlo' distribution"); rug(xs)

x <- seq(from = -10, to = 1010, by = 2)
plot(x = x, y = dmaxlo(x, shape = 4, scale = 1000),
     type = "l", ylab = "dens",
     col = "orangered", main = "dmaxlo and dgpd")
abline(h = 0)
lines(x = x, y = dgpd(x, shape = -1/4, scale = 250),
     type = "l",
     col = "SpringGreen3", lty = "dashed")



}
\keyword{distribution}

