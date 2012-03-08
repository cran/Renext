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
  distribution. 
}
\usage{
   dlomax(x, scale = 1.0, shape = 4.0, log = FALSE)
   plomax(q, scale = 1.0, shape = 4.0, lower.tail = FALSE)
   qlomax(p, scale = 1.0, shape = 4.0)
   rlomax(n, scale = 1.0, shape = 4.0) 
}
\arguments{
  \item{x, q}{Vector of quantiles.}
  \item{p}{Vector of probabilities.}
  \item{n}{Number of observations.}
  \item{scale, shape}{Shift and shape parameters. Vectors of length
    > 1 are not accepted.}
  \item{log}{Logical; if \code{TRUE}, the log density is returned.}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities
    are \eqn{\textrm{Pr}[X <= x]}{Pr[X <= x]}, otherwise,
    \eqn{\textrm{Pr}[X > x]}{X>x}}.
}
\details{
  The lomax distribution function with shape \eqn{\rho>0}{rho>0} and
  scale \eqn{\delta>0}{delta>0} has distribution function
  \deqn{F(y) = 1- \left[\frac{\delta}{\delta+y}\right]^\rho \quad y>0}{
    F(y) = 1 - [delta /(delta+y)]^rho      y > 0}
  This distribution has increasing hazard and decreasing mean residual
  life (MRL). The coefficient of variation decreases with
  \eqn{\rho}{rho} and tends to \eqn{1}{1} for large \eqn{\rho}{rho}.
  The default value \eqn{\rho=4}{rho = 4} corresponds to \eqn{CV = \sqrt{2}}{CV
    = sqrt(2)}.
}
\value{
  \code{dlomax} gives the density function, \code{plomax} gives the
  distribution function, \code{qlomax} gives the quantile function,
  and \code{rlomax} generates random deviates.
}
\note{
  This distribution is sometimes called \emph{log-exponential}
  and also called the \emph{Lomax distribution}. It
  is occasionally used in POT with \code{scale} taken as
  the threshold. 
}
\references{
  Johnson N. Kotz S. and N. Balakrishnan 
  \emph{Continuous Univariate Distributions} vol. 1, Wiley 1994.
}

%%\references{  }	
%%\seealso{\code{\link{}} }

\examples{
library(Renext)
shape <- 5; scale <- 10

xl <- qlomax(c(0.00, 0.99), scale = scale, shape = shape)
x <- seq(from = xl[1], to = xl[2], length.out = 200)
f <- dlomax(x, scale = scale, shape = shape)
plot(x, f, type = "l", main = "Lomax distribution density")

F <- plomax(x, scale = scale, shape = shape)
plot(x, F, type ="l", main ="Lomax distribution function")
}
\keyword{distribution}

