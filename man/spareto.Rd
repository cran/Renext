\name{spareto}
\alias{spareto}
\alias{dspareto}
\alias{pspareto}
\alias{qspareto}
\alias{rspareto}
\title{Shifted Pareto dstribution}
\description{
  Density function, distribution function, quantile function and
  random generation for the shifted Pareto distribution. 
  distribution. 
}
\usage{
   dspareto(x, delta = 1.0, shape = 4.0, log = FALSE)
   pspareto(q, delta = 1.0, shape = 4.0, lower.tail = FALSE)
   qspareto(p, delta = 1.0, shape = 4.0)
   rspareto(n, delta = 1.0, shape = 4.0) 
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{delta, shape}{Shift and shape parameters. Vectors of length
    > 1 are not accepted.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities
    are \eqn{\textrm{Pr}[X <= x]}{Pr[X <= x]}, otherwise,
    \eqn{\textrm{Pr}[X > x]}{X>x}}.
}
\details{
  The spareto distribution function with shape \eqn{\rho>0}{rho>0} and
  shift \eqn{\delta>0}{delta>0} has distribution function
  \deqn{F(y) = 1- \left[\frac{\delta}{\delta+y}\right]^\rho \quad y>0}{
    F(y) = 1 - [delta /(delta+y)]^rho      y > 0}
  This distribution has increasing hazard and decreasing mean residual
  life (MRL). The coefficient of variation decreases with
  \eqn{\rho}{rho} and tends to \eqn{1}{1} for large \eqn{\rho}{rho}.
  The default value \eqn{\rho=4}{rho = 4} corresponds to \eqn{CV = \sqrt{2}}{CV
    = sqrt(2)}.
  
  This distribution is sometimes called \emph{log-exponential}. It
  is occasionally used in POT with the shift \code{delta} taken as
  the threshold. 
}
\value{
  \code{dspareto} gives the density function, \code{pspareto} gives the
  distribution function, \code{qspareto} gives the quantile function,
  and \code{rspareto} generates random deviates.
}
%%\references{  }	
%%\seealso{\code{\link{}} }

\examples{
library(Renext)
shape <- 5; delta <- 10

xl <- qspareto(c(0.00, 0.99), delta = delta, shape = shape)
x <- seq(from = xl[1], to = xl[2], length.out = 200)
f <- dspareto(x, delta = delta, shape = shape)
plot(x, f, type = "l", main ="Shifted Pareto distribution density")

F <- pspareto(x, delta = delta, shape = shape)
plot(x, F, type ="l", main ="Shifted Pareto distribution function")
}
\keyword{distribution}

