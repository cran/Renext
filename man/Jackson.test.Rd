\name{Jackson.test}
\alias{Jackson.test}
\title{Jackson's test of exponentiality}
\usage{
   Jackson.test(x, method = c("num", "sim", "asymp"), nSamp = 15000)
}
\arguments{
  \item{x}{numeric vector or matrix.}
  
  \item{method}{Character: choice of the method used to
    compute the \eqn{p}-value. See the \bold{Details}
    section.}
  
  \item{nSamp}{Number of samples used to compute the
    \eqn{p}-value if \code{method} is \code{"sim"}.}
}
\value{

  A list of results.
  
  \item{statistic, p.value}{

    The statistic and \eqn{p}-value.

  }
  \item{df}{

    Number \eqn{n} of observations.

  }
  \item{method}{

    Description of the test implemented, regardless of how the
    \eqn{p}-value has been computed.

  } 
 
  
}
\description{
  
  Jackson's test of exponentiality

}
\details{
  
  Compute the Jackson's test of exponentiality.  The test
  statistic is the ratio of weighted sums of the order
  statistics. Both sums can also be written as weighted sums
  of the \emph{scalings}.
  
  The Jackson's statistic for a sample of size \eqn{n} of the
  exponential distribution can be shown to be approximately
  normal. More precisely \eqn{\sqrt{n}(J_n -2)}{sqrt(n)*(Jn
    -2)} has approximately a standard normal distribution.
  This distribution is used to compute the \eqn{p}-value when
  \code{method} is \code{"asymp"}. When \code{method} is
  \code{"num"}, a numerical approximation of the distribution
  is used. Finally, when \code{method} is \code{"sim"} the
  \eqn{p}-value is computed by simulating \code{nSamp}
  samples of size \code{length(x)} and estimating the
  probability to have a Jackson's statistic larger than that
  of the 'observed' \code{x}.

}
\note{
  
  Jackson's test of exponentiality works fine for a Lomax alternative
  (GPD with heavy tail). It then reaches nearly the same power as a
  Likelihood Ratio (LR) test, see Kozubowski et al.  It can be
  implemented more easily than the LR test because simulated values of
  the test statistic can be obtained quickly enough to compute the
  \eqn{p}-value by simulation.
  
}
\author{

  Yves Deville

}
\seealso{

  The \code{\link{Jackson}} function computing the statistic and the
  \code{\link{LRExp.test}} function.

}
\references{
  
  J. Beirlant and T. de Weit and Y. Goegebeur(2006) "A Goodness-of-fit
  Statistic for Pareto-Type Behaviour", \emph{J. Comp. Appl. Math.},
  186(1), pp. 99-116.

  T.J. Kozubowski, A. K. Panorska, F. Qeadan, A. Gershunov and
  D. Rominger (2009)
  "Testing Exponentiality Versus Pareto Distribution via Likelihood Ratio"
  \emph{Comm. Statist. Simulation Comput.} 38(1), pp. 118-139.
  
}
\examples{
set.seed(1234)
x <- rGPD(n = 50, loc = 0, scale = 1, shape = 0.1)
Jackson.test(x, method = "num")$p.value
Jackson.test(x, method = "asymp")$p.value
Jackson.test(x, method = "sim")$p.value
}
